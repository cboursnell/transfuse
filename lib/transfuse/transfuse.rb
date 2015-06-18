class String
  def revcomp
    self.tr("ACGT", "TGCA").reverse
  end
end

module Transfuse

  require 'csv'
  require 'transrate'
  require 'threach'

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

    def sequence_alignment clusters
      output_files = {}
      clusters.threach(@threads) do |id, list| # threach
        output_file = "consensus_output_#{Thread.current.object_id}.fa"
        output_files[output_file]=1
        File.open(output_file, "ab") do |out|
          # puts "msa: #{id}\tsequences: #{list.length}" if id == "48"
          if list.size > 1
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
            # File.open("cluster#{id}.fa", "wb") { |out| out.write seq }
            output = "cluster#{id}.aln.fa"
            print "Thread #{Thread.current.object_id} clustalo on cluster #{id}...\n" if @verbose
            cmd = "echo  \"#{seq}\" | #{@clustalo} -i - -o #{output} "
            cmd << "--outfmt clu --force "
            cmd << "--output-order tree-order --infmt fa --wrap 50000"
            align = Cmd.new cmd
            align.run

            unless align.status.success?
              abort align.stderr
            end
            consensus = parse_msa(output)
            out.write ">cluster#{id}\n"
            out.write "#{consensus}\n"
            File.delete(output)
          else
            out.write ">cluster#{id}\n"
            out.write "#{@sequences[list[0][:name]]}\n"
          end
        end
      end
      cmd = "cat "
      cmd << output_files.join(" ")
      cmd << " > consensus_output.fa"
      cat = Cmd.new cmd
      cat.run
    end

    def parse_msa(output)
      msa = {}
      exons = {}
      longest = 0
      File.open(output).each do |line|
        if line =~ /CLUSTAL/ or line.length < 2 or line=~/^\s+/
        else
          if line =~ /(.+)\s+(.+)/
            name = $1
            alignment = $2
            msa[name] = []
            prev = "-"
            e = 0
            alignment.each_char.with_index do |c, index|
              msa[name][index] = c
              longest = index if index > longest
              if prev == "-" and c != "-"
                e += 1
              end
              prev = c
            end
            # puts "#{name}\t#{e}"
            exons[name] = e
          end
        end
      end
      points=[]
      column = ""
      consensus = ""
      (0..longest).each do |i|
        column = ""
        msa.each do |name, list|
          if exons[name]==1
            if list[i] != "-"
              column << list[i]
            end
          end
        end
        highest=-1
        best=0
        bestbase="N"
        column.split("").uniq.each_with_index do |base, index|
          if column.count(base) > highest
            best = index
            bestbase = base
          end
        end
        if bestbase == "N"
          puts "output: #{output}"
          puts "column: #{column}\t#{i}"
          p exons
          exit
        end
        consensus << bestbase

      end
      return consensus
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
        if !File.exist?(new_filename) or File.stat(new_filename).size < 1
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
