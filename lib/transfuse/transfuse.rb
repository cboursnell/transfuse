class String
  def revcomp
    self.tr("ACGT", "TGCA").reverse
  end
end

module Transfuse

  require 'bio'
  require 'csv'
  require 'transrate'
  require 'threach'

  class Transfuse

    def initialize threads, verbose
      @threads = threads
      @verbose = verbose
    end

    def check_dependencies
      # Check dependencies if they are relevant to the command issued,
      # and handle any commands to install missing ones
      gem_dir = Gem.loaded_specs['transfuse'].full_gem_path
      gem_deps = File.join(gem_dir, 'deps', 'deps.yaml')

      return Bindeps.missing gem_deps

    end # check_dependencies

    def install_dependencies
      # Check dependencies if they are relevant to the command issued,
      # and handle any commands to install missing ones
      gem_dir = Gem.loaded_specs['transfuse'].full_gem_path
      gem_deps = File.join(gem_dir, 'deps', 'deps.yaml')

      Bindeps.require gem_deps

    end # check_dependencies

    def check_files string, option
      # puts "check file string: #{string}" if @verbose
      abort "Please specify --#{option} option" if string.nil?
      list = []
      string.split(",").each do |file|
        file = File.expand_path(file)
        if File.exist?(file)
          puts "#{File.basename(file)} exists" if @verbose
          list << file
        else
          abort "#{File.basename(file)} not found"
        end
      end
      return list
    end

    def concatenate assemblies
      catted_fasta = "all-"
      fasta = []
      assemblies.each do |name|
        fasta << File.basename(name, File.extname(name))[0..5]
      end
      catted_fasta << fasta.join("-")
      catted_fasta << ".fa"
      puts "concatenating assemblies into #{catted_fasta}" if @verbose
      cmd = "cat "
      assemblies.each do |file|
        cmd << " #{file} "
      end
      cmd << " > #{catted_fasta}"
      catter = Cmd.new cmd
      catter.run catted_fasta
      return File.expand_path(catted_fasta)
    end

    def load_fasta fasta
      print "loading fasta sequence #{fasta}..." if @verbose
      @sequences = {}
      count = 1
      Bio::FastaFormat.open(fasta).each do |entry|
        @sequences[entry.entry_id] = entry.seq.to_s
        print "." if count%10_000==0 and @verbose
        count +=1
      end
      puts " Done" if @verbose
    end

    def cluster file, id
      puts "clustering #{file}" if @verbose
      cluster = Cluster.new @threads, @verbose, id
      return cluster.run file
    end

    def consensus msa, scores, output
      cons = Consensus.new(@verbose)
      return cons.run(msa, scores, output)
    end

    def load_scores files
      scores = {}
      files.each do |file|
        CSV.foreach(file, :headers => true,
                          :header_converters => :symbol,
                          :converters => :all) do |row|
          name = row[:contig_name]
          scores[name] = { :score => row[:score].to_f,
                           :p_good => row[:p_good].to_f,
                           :p_bases_covered => row[:p_bases_covered].to_f,
                           :coverage => row[:coverage].to_f }
        end
      end
      return scores
    end

    def filter files, scores
      filtered_files = []
      files.each_with_index do |file, index|
        new_filename = "#{File.basename(file, File.extname(file))}_filtered.fa"
        if !File.exist?(new_filename) or File.stat(new_filename).size < 1
          File.open(new_filename, "wb") do |out|
            puts "filtering #{file}..." if @verbose
            Bio::FastaFormat.open(file).each do |entry|
              contig_name = entry.entry_id
              contig_name = "contig#{index}_#{contig_name}"
              if scores.key?(contig_name) and
                 scores[contig_name][:score] > 0.01 and
                 scores[contig_name][:coverage] >= 1
                out.write ">#{contig_name}\n"
                out.write "#{entry.seq}\n"
              elsif !scores.key?(contig_name)
                abort "Can't find '#{contig_name}' in scores"
              end
            end
          end
        end
        filtered_files << File.expand_path(new_filename)
      end
      return filtered_files
    end

    def transrate_consensus file, output, left, right
      output = File.expand_path(output)
      puts "transrate on #{file}" if @verbose
      file = File.expand_path(file)
      name = File.basename(file, File.extname(file))
      dir = "transrate_#{name}"
      Dir.mkdir(dir) unless Dir.exist?(dir)
      Dir.chdir(dir) do
        assembly = Transrate::Assembly.new(file)
        transrater = Transrate::Transrater.new(assembly, nil, threads:@threads)
        rename = "assembly_#{name}_score_optimisation.csv"
        rm = transrater.read_metrics(left.join(','), right.join(','))
        stats = rm.read_stats
        File.rename("assembly_score_optimisation.csv", rename)
        scores={}
        assembly.each do |name, contig|
          scores[name] = { :score => contig.score.to_f,
                           :p_good => contig.p_good.to_f,
                           :p_bases_covered => contig.p_bases_covered.to_f,
                           :coverage => contig.coverage.to_f }
        end
        scores_file = "#{name}_scores.csv"
        stats_file = "../#{name}_stats.txt"
        puts "  writing scores" if @verbose
        File.open(scores_file, "wb") do |out|
          scores.each do |name, hash|
            out.write "#{name}\t#{hash[:score]}\t#{hash[:p_good]}\t"
            out.write "#{hash[:p_bases_covered]}\t#{hash[:coverage]}\n"
          end
        end
        puts "  writing filtered fasta file" if @verbose
        File.open(output, "wb") do |out|
          assembly.each do |name, contig|
            if contig.score.to_f > 0.01 and contig.coverage.to_f >= 1
              out.write ">#{name}\n"
              out.write "#{contig.seq.seq}\n"
            end
          end
        end
        puts "  writing stats" if @verbose
        File.open(stats_file, "wb") do |out|
          stats.each do |key, value|
            out.write "#{key}\t#{value}\n"
          end
          out.write "assembly score:\t#{transrater.assembly_score}\n"
          optimal = transrater.assembly_optimal_score("prefix")
          out.write "optimal score :\t#{optimal[0]}\n"
          out.write "cutoff        :\t#{optimal[1]}\n"
        end
      end
    end

    def transrate files, left, right
      unless left.is_a?(Array)
        left = [left]
      end
      unless right.is_a?(Array)
        right = [right]
      end
      scores = {}
      shortname = ""
      files.each do |n|
        a = File.basename(n).split("_").first
        if a.length > 5
          a = a[0..4]
        end
        # shortname << File.basename(n, File.extname(n))[0..4]
      end
      scores_file = "#{shortname}_scores.csv"
      if File.exist?(scores_file)
        puts "loading scores from file" if @verbose
        File.open(scores_file).each do |line|
          name, score, p_good, p_bases_covered, coverage = line.chomp.split("\t")
          scores[name] = { :score => score.to_f,
                           :p_good => p_good.to_f,
                           :p_bases_covered => p_bases_covered.to_f,
                           :coverage => coverage.to_f }
        end
      else
        files.each_with_index do |fasta, index|
          puts "transrate on #{fasta}" if @verbose
          dir = "transrate_#{File.basename(fasta, File.extname(fasta))}"
          Dir.mkdir(dir) unless Dir.exist?(dir)
          Dir.chdir(dir) do
            assembly = Transrate::Assembly.new(fasta)
            transrater = Transrate::Transrater.new(assembly, nil, threads:@threads)
            rename = "assembly#{index}_score_optimisation.csv"
            transrater.read_metrics(left.join(','), right.join(','))
            File.rename("assembly_score_optimisation.csv", rename)
            assembly.each do |name, contig|
              name = "contig#{index}_#{name}"
              scores[name] = { :score => contig.score.to_f,
                               :p_good => contig.p_good.to_f,
                               :p_bases_covered => contig.p_bases_covered.to_f,
                               :coverage => contig.coverage.to_f }

            end
            File.open("summary.txt","w") do |out|
              out.write "fasta\tscore\toptimal\tcutoff\n"
              out.write "#{fasta}\t#{transrater.assembly_score}\t#{transrater.assembly_optimal_score("prefix").join("\t")}\n"
            end
          end
        end
        File.open(scores_file, "wb") do |out|
          scores.each do |name, hash|
            out.write "#{name}\t#{hash[:score]}\t#{hash[:p_good]}\t"
            out.write "#{hash[:p_bases_covered]}\t#{hash[:coverage]}\n"
          end
        end
      end
      return scores
    end

  end

end
