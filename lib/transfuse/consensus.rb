
require 'bio'
require 'set'

module Transfuse

  class Consensus

    attr_reader :contigs

    def initialize verbose
      @verbose = verbose
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
          msa.each do |id, list|
            count+=1
            print "." if count%5_000==0 and @verbose
            exons={}
            cons = {}
            length = list[0][:seq].length
            list.each_with_index do |hash, index|
              seq = hash[:seq]
              name = hash[:name]
              out2.write "#{id}\t#{scores[name][:score]}\t#{name}\n"
              prev = ""
              # gap = 0
              exon = 0
              seq.each_char do |c|
                if c=="-"
                  base="-"
                else
                  base="*"
                end
                if base!=prev
                  if c=="-"
                    # gap+=1
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
              exons[index] = exon
            end

            consensus = ""
            0.upto(length-1) do |i|
              base = "N"
              counts = {}
              list.each_with_index do |hash, index|
                seq = hash[:seq]
                if seq[i] != "-" and seq[i] != "N"
                  counts[seq[i]] ||= 0
                  counts[seq[i]] += 1
                  if exons[index] == 1
                    base = seq[i]
                  end
                end
              end
              if counts.size>0
                base = counts.sort.last.first
              end
              consensus << base
            end

            if consensus.count("N") < consensus.length.to_f*0.5
              cons[consensus] = 1
            end

            list.each_with_index do |hash, index|
              if exons[index] > 1
                seqn = hash[:seq].delete("-")
                cons[hash[:seq].delete("-")] = 1
              end
            end

            cons.each_with_index do |s,index|
              out.write ">contig#{id}.#{index+1}\n"
              out.write "#{s[0]}\n"
            end

          end # msa.each
        end # file
      end # file open
      puts " Done" if @verbose
      return preoutput
    end # def

  end

end