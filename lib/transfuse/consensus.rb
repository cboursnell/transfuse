
require 'bio'
require 'set'

module Transfuse

  class Consensus

    attr_reader :contigs

    def initialize verbose
      @verbose = verbose
    end

    def count_exons list
      exon_counts = {}
      list.each_with_index do |hash, index|
        seq = hash[:seq]
        exon_count=0
        gap_count=0
        prev=""
        seq.each_char do |c|
          if c=="-"
            base = "-"
          else
            base = "*"
          end

          if base!=prev
            if c=="-"
              gap_count+=1
            else
              exon_count+=1
            end
          end

          if c=="-"
            prev = "-"
          else
            prev = "*"
          end
        end
        exon_counts[exon_count] ||= []
        exon_counts[exon_count] << index
      end
      return exon_counts
    end

    def run msa, scores, output
      preoutput = "#{File.basename(output, File.extname(output))}_cons.fa"
      return preoutput if File.exist?(preoutput)
      print "writing consensus " if @verbose
      # msa is a hash
      #   key = cluster id
      #   value = list
      #     list of sequences in cluster aligned with gaps
      count = 0
      File.open("#{output}.data", "w") do |out2|
        File.open(preoutput, "w") do |out|
          msa.each do |id, seq_list|
            count+=1
            print "." if count%5_000==0 and @verbose
            exons={}
            cons = {}
            length = seq_list[0][:seq].length
            exons = count_exons(seq_list)

            exons.each do |count, list|
              0.upto(length-1) do |pos|
                base = "N"
                list.each do |index|
                  b = seq_list[index][:seq][pos]
                  if b != "-" and b != "N"
                    base = b
                  end
                end
                if base != "N"
                  cons[count]||=""
                  cons[count]<<base
                end
              end
            end

            cons.each_with_index do |s, index|
              out.write ">contig#{id}.#{index+1}\n"
              out.write "#{s[1]}\n"
            end

          end # msa.each
        end # file
      end # file open
      puts " Done" if @verbose
      return preoutput
    end # def

  end

end