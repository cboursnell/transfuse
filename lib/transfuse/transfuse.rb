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

    def load_fasta fasta
      puts "loading fasta sequence #{fasta}" if @verbose
      @sequences = {}
      Bio::FastaFormat.open(fasta).each do |entry|
        @sequences[entry.entry_id] = entry.seq.to_s
      end
    end

    def cluster file
      puts "clustering #{file}" if @verbose
      cluster = Cluster.new @threads, @verbose
      return cluster.run file
    end

    def consensus clusters
      # TODO add threach to this function
      clusters.each do |id, list|
        puts "#{id}"
        con = Consensus.new 31
        # TODO make this check that the cluster has more than one sequence in
        #      it before running consensus
        list.each_with_index do |hash,index|
          con.add_kmers index, @sequences[hash[:name]]
        end
        con.output id # TODO make this get output and write to files
      end
      # TODO concatenate all the individual fasta sequences
      #      with one giant 'cat' command?
      seqs
    end

    def load_scores files
      scores = {}
      files.each do |file|
        CSV.foreach(file, :headers => true,
                          :header_converters => :symbol,
                          :converters => :all) do |row|
          name = row[:contig_name]
          scores[name] = { :score => row[:score],
                           :p_good => row[:p_good],
                           :p_bases_covered => row[:p_bases_covered],
                           :coverage => row[:coverage] }
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
            puts "opening #{file}..."
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
          File.rename("assembly_score_optimisation.csv", "assembly#{index}_score_optimisation.csv")
          assembly.each do |name, contig|
            name = "contig#{index}_#{name}"
            scores[name] = { :score => contig.score,
                             :p_good => contig.p_good,
                             :p_bases_covered => contig.p_bases_covered,
                             :coverage => contig.coverage }

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
