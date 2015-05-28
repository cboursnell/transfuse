class String
  def revcomp
    self.tr("ACGT", "TGCA").reverse
  end
end

module Transfuse

  require 'csv'
  require 'transrate'

  class Transfuse

    def initialize threads, verbose
      @threads = threads
      @verbose = verbose
      @clustalo = Which::which('clustalo').first
      raise "clustalo was not in the PATH - please install it" unless @clustalo
    end

    def check_files string
      list = []
      string.split(",").each do |file|
        file = File.expand_path(file)
        if File.exist?(file)
          puts "#{file} exists" if @verbose
          list << file
        else
          abort "#{file} not found"
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

    def cluster file
      puts "clustering #{file}" if @verbose
      cluster = Cluster.new @threads, @verbose
      return cluster.run file
    end

    def load_fasta fasta
      @sequences = {}
      Bio::FastaFormat.open(fasta).each do |entry|
        @sequences[entry.entry_id] = entry.seq.to_s
      end
    end

    def sequence_alignment clusters
      clusters.each do |id, list| # threach
        if list.size > 5
          seq = ""
          list.each do |hash|
            seq << ">#{hash[:name]}\n"
            if hash[:strand] == "+"
              seq << "#{@sequences[hash[:name]]}\n"
            elsif hash[:strand] == "-"
              seq << "#{@sequences[hash[:name]].revcomp}\n"
            else
              abort "Unknown strand #{hash[:strand]}"
            end
          end
          cmd = "echo -e \"#{seq}\" | #{@clustalo} -i - --outfmt fa "
          cmd << "--output-order tree-order"
          align = Cmd.new cmd
          align.run
          File.open("cluster#{id}.fa", "wb") do |out|
            out.write align.stdout
          end
        end
      end
    end

    def load_scores files
      scores = {}
      files.each do |file|
        CSV.foreach(file, :headers => true,
                          :header_converters => :symbol,
                          :converters => :all) do |row|
          name = row[:contig_name]
          score = row[:score]
          scores[name] = score
        end
      end
      return scores
    end

    def filter files, scores
      filtered_files = []
      files.each_with_index do |file, index|
        new_filename = "#{File.basename(file, File.extname(file))}_filtered.fa"
        unless File.exist?(new_filename)
          File.open(new_filename, "wb") do |out|
            puts "opening #{file}..."
            Bio::FastaFormat.open(file).each do |entry|
              contig_name = entry.entry_id
              contig_name = "contig#{index}_#{contig_name}"
              if scores.key?(contig_name) and scores[contig_name] > 0.01
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


    def transrate files, left, right
      scores = {}
      scores_file = "scores.csv"
      if File.exist?(scores_file)
        puts "loading scores from file" if @verbose
        File.open(scores_file).each do |line|
          name, score = line.chomp.split("\t")
          scores[name] = score.to_f
        end
      else
        files.each_with_index do |fasta, index|
          puts "transrate on #{fasta}" if @verbose
          assembly = Transrate::Assembly.new(fasta)
          transrater = Transrate::Transrater.new(assembly, nil, threads:@threads)
          transrater.read_metrics(left.join(','), right.join(','))
          assembly.each do |name, contig|
            name = "contig#{index}_#{name}"
            scores[name] = contig.score
          end
        end
        File.open(scores_file, "wb") do |out|
          scores.each do |name, score|
            out.write "#{name}\t#{score}\n"
          end
        end
      end
      return scores
    end

    def select_contigs clusters, scores
      puts "selecting contigs" if @verbose
      best = []
      clusters.each do |cluster_id, list|
        best_score = 0
        best_contig = ""
        list.each do |contig_name|
          unless scores[contig_name]
            abort "can't find #{contig_name} in scores hash\n"
          end
          if scores[contig_name] > best_score
            best_score = scores[contig_name]
            best_contig = contig_name
          end
        end
        best << best_contig
      end
      return best
    end

    def output_contigs best, fasta, output
      puts "writing contigs" if @verbose
      # read in catted fasta sequences
      sequences = {}
      Bio::FastaFormat.open(fasta).each do |entry|
        sequences[entry.entry_id] = entry.seq
      end
      File.open(output, "wb") do |out|
        best.each do |contig_name|
          if sequences.key?(contig_name)
            out.write ">#{contig_name}\n"
            out.write "#{sequences[contig_name]}\n"
          else
            puts "can't find #{contig_name} in #{fasta}"
          end
        end
      end
    end

  end

end
