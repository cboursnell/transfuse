module Transfuse

  require 'csv'

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
      output = cluster.run file
      return cluster.parse_output(output)
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

    def transrate files, left, right
      rate = Transrate.new @threads
      score_files = rate.run files, left, right
      return load_scores(score_files)
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
