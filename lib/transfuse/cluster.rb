module Transfuse

  require 'bio'
  require 'fixwhich'

  class Cluster

    def initialize threads
      @vsearch = Which::which('vsearch').first
      raise RuntimeError unless @vsearch
      @id = "1.00"
      @threads = threads
    end

    def run fasta
      vsearch fasta
    end

    def vsearch fasta
      cluster_output = "clusters.txt"
      vsearch_cmd = generate_command fasta, cluster_output
      cluster = Cmd.new vsearch_cmd
      cluster.run cluster_output
      return cluster_output
    end

    def generate_command fasta, out
      vsearch = "#{@vsearch}"
      vsearch << " --cluster_fast #{fasta}"
      vsearch << " --id #{@id}"
      vsearch << " --strand both"
      vsearch << " --uc #{out}"
      vsearch << " --threads #{@threads}"
      return vsearch
    end

    def parse_output cluster_output
      clusters = {}
      File.open(cluster_output).each_line do |line|
        if line.start_with?("S") or line.start_with?("H")
          cols = line.chomp.split("\t")
          cluster = cols[1].to_i
          contig_name = cols[8]
          clusters[cluster] ||= []
          clusters[cluster] << contig_name
        end
      end
      return clusters
    end

  end

end