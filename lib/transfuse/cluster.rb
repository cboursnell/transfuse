module Transfuse

  require 'bio'
  require 'fixwhich'

  class Cluster

    def initialize threads, verbose, id
      @vsearch = Which::which('vsearch').first
      raise "vsearch was not in the PATH - please install it" unless @vsearch
      @id = id.to_s
      @threads = threads
      @verbose = verbose
    end

    def run fasta
      cluster_output, msa_output = vsearch fasta
      return parse_vsearch_output(cluster_output, msa_output)
    end

    def vsearch fasta
      print "running vsearch..." if @verbose
      cluster_output = "#{File.basename(fasta)}-#{@id}.clust"
      msa_output = "#{File.basename(fasta)}-#{@id}.aln"
      vsearch_cmd = generate_vsearch_command fasta, cluster_output, msa_output
      cluster = Cmd.new vsearch_cmd
      cluster.run cluster_output
      puts " Done. Created #{cluster_output}" if @verbose
      return [cluster_output, msa_output]
    end

    def generate_vsearch_command fasta, out, msa
      vsearch = "#{@vsearch}"
      vsearch << " --cluster_fast #{fasta}"
      vsearch << " --id #{@id}"
      vsearch << " --iddef 0" # cd-hit definition of sequence id
      vsearch << " --qmask none" # no masking
      vsearch << " --strand both"
      vsearch << " --uc #{out}"
      vsearch << " --msaout #{msa}"
      vsearch << " --threads #{@threads}"
      return vsearch
    end

    def parse_vsearch_output cluster_output, msa_output
      print "parsing vsearch output" if @verbose
      clusters = {}
      lookup = {}
      second = 0
      count = 0
      File.open(cluster_output).each_line do |line|
        count+=1
        if line.start_with?("S") or line.start_with?("H")
          cols = line.chomp.split("\t")
          cluster = cols[1]
          len = cols[2].to_i
          cigar = cols[7]
          strand = cols[4]
          strand = "+" if strand == "*"
          contig_name = cols[8]

          clusters[cluster] ||= []
          clusters[cluster] << { :name => contig_name, :strand => strand }
          lookup[contig_name] = cluster
        end
        if count%10_000==0 and @verbose
          print "."
        end
      end
      puts " Done" if @verbose
      print "parsing msa output    " if @verbose
      count = 0
      msa = {}
      Bio::FastaFormat.open(msa_output).each do |entry|
        count += 1
        name = entry.entry_id
        if name != "consensus"
          # name = name[1..-1]
          if name[0]=="*"
            name = name[1..-1]
          end
          # what cluster is name in?
          cluster = lookup[name]
          msa[cluster] ||= []
          msa[cluster] << { :name => name, :seq => entry.seq.seq }
        end
        if count%10_000==0 and @verbose
          print "."
        end

      end
      puts " Done" if @verbose
      return msa
    end

  end

end
