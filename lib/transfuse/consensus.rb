
require 'bio'
require 'set'

module Transfuse

  class Consensus

    attr_reader :contigs

    def initialize
    end

    def run msa, output
      # msa is a hash
      #   key = cluster id
      #   value = list
      #      list of sequences in cluster aligned with gaps
      File.open(output, "w") do |out|
        msa.each do |id, list|
          exons={}
          cons = []
          length = list[0].length
          list.each_with_index do |seq, index|
            prev = ""
            gap=0
            exon=0
            seq.each_char do |c|
              if c=="-"
                base="-"
              else
                base="*"
              end
              if base!=prev
                if c=="-"
                  gap+=1
                else
                  exon+=1
                end
              end
              if c=="-"
                prev = "-"
              else
                prev = "*"
              end
            end
            exons[index]=exon
          end

          consensus = ""
          0.upto(length-1) do |i|
            base="N"
            list.each_with_index do |seq, index|
              if exons[index]==1
                if seq[i] != "-"
                  base = seq[i]
                end
              end
            end
            consensus << base
          end

          cons << consensus

          list.each_with_index do |seq, index|
            if exons[index]>1
              cons << seq.delete("-")
            end
          end

          cons.each_with_index do |s,index|
            out.write ">contig#{id}.#{index+1}\n"
            out.write "#{s}\n"
          end

        end # msa.each
      end # file
    end # def

  end

end