module Transfuse

  require 'csv'
  require 'transrate'

  class Transfuse

    def initialize threads
      @threads = threads
    end

    def check_files string
      list = []
      string.split(",").each do |file|
        file = File.expand_path(file)
        if File.exist?(file)
          list << file
        else
          abort "#{file} not found"
        end
      end
      return list
    end

    def concatenate assemblies
      catted_fasta = "all.fa"
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
      cluster = Cluster.new @threads
      return cluster.run file
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
        File.open(new_filename, "wb") do |out|
          Bio::FastaFormat.open(file).each do |entry|
            contig_name = entry.entry_id
            if scores.key?(contig_name) and scores[contig_name] > 0.01
              out.write ">#{index}_#{contig_name}\n"
              out.write "#{entry.seq}\n"
            end
          end
        end
        filtered_files << File.expand_path(new_filename)
      end
      return filtered_files
    end

    def transrate files, left, right
      scores = {}
      files.each do |fasta|
        assembly = Transrate::Assembly.new(fasta)
        transrater = Transrate::Transrater.new(assembly, nil, threads:@threads)
        transrater.read_metrics(left.join(','), right.join(','))
        assembly.each do |name, contig|
          scores[name] = contig.score
        end
      end
      return scores
    end

    def select_contigs clusters, scores
      best = []
      clusters.each do |cluster_id, list|
        best_score = 0
        best_contig = ""
        list.each do |contig_name|
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
