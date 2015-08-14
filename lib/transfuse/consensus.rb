
require 'bio'
require 'set'

module Transfuse

  class Consensus

    attr_reader :contigs

    def initialize kmer
      @kmer = kmer # size of kmers
      @contigs = {}
      @match = 2
      @mismatch = -10
      @gap_open = 3
      @gap_extend = 1
    end

    def add_fasta file
      puts "opening #{file}"
      id = 0
      Bio::FastaFormat.open(file).each do |entry|
        add_kmers id, entry.seq.seq
        id += 1
      end
    end

    def add_kmers id, seq
      types = {}
      cons = {}
      if id > 0
        0.upto(id-1) do |c|
          unless @contigs[c].nil?
            type, consensus = compare(kmerise(c, @contigs[c]),
                                      kmerise(id, seq) )
            types[c] = type
            cons[c] = consensus

            if type > 0 and type <= 4
              @contigs.delete(c)
              seq = consensus
            end
          end
        end
      end
      @contigs[id]=seq
    end

    def output id
      seq = ""
      @contigs.each do |k, seq|
        seq << ">contig#{id}.#{k}\n"
        seq << "#{seq}\n"
      end
      return seq
    end

    def compare kmers1, kmers2
      s1 = kmers1.length
      s2 = kmers2.length
      matrix = initialize_matrix s1, s2

      1.upto(s1) do |x|
        1.upto(s2) do |y|
          next if x==0 and y==0

          # m matrix
          if x > 0 and y > 0
            if kmers1[x-1] == kmers2[y-1] # match
              score_m, arrow_m = max3(matrix[x-1][y-1][:score_m] + @match,
                                    matrix[x-1][y-1][:score_x] + @match,
                                    matrix[x-1][y-1][:score_y] + @match)

            else
              score_m, arrow_m = max3(matrix[x-1][y-1][:score_m] + @mismatch,
                                    matrix[x-1][y-1][:score_x] + @mismatch,
                                    matrix[x-1][y-1][:score_y] + @mismatch)
            end
          end

          # x matrix
          if x > 0
            score_x, arrow_x = max3(matrix[x-1][y][:score_m] - @gap_open,
                                    matrix[x-1][y][:score_x] - @gap_extend,
                                    matrix[x-1][y][:score_y] - @gap_open )
          end

          # y matrix
          if y > 0
            score_y, arrow_y = max3(matrix[x][y-1][:score_m] - @gap_open,
                                    matrix[x][y-1][:score_x] - @gap_open,
                                    matrix[x][y-1][:score_y] - @gap_extend )
          end

          matrix[x][y] = { :score_m => score_m,
                           :score_x => score_x,
                           :score_y => score_y,
                           :arrow_m => arrow_m,
                           :arrow_x => arrow_x,
                           :arrow_y => arrow_y }
        end # y
      end # x

      x = s1
      y = s2

      alignment1=[]
      alignment2=[]

      score, layer = max3(matrix[x][y][:score_m],
                          matrix[x][y][:score_x],
                          matrix[x][y][:score_y])
      destination_layer = ""

      while x > 0 or y > 0
        if layer =~ /M/
          destination_layer = matrix[x][y][:arrow_m]
          alignment1 << kmers1[x-1]
          alignment2 << kmers2[y-1]
          x-=1
          y-=1
        elsif layer =~ /X/
          destination_layer = matrix[x][y][:arrow_x]
          alignment1 << kmers1[x-1]
          alignment2 << "-"
          x-=1
        elsif layer =~ /Y/
          destination_layer = matrix[x][y][:arrow_y]
          alignment1 << "-"
          alignment2 << kmers2[y-1]
          y-=1
        else

        end
        layer = destination_layer

      end # while

      a = {:gaps => 0, :exons => 0, :overlaps => 0}
      b = {:gaps => 0, :exons => 0, :overlaps => 0}

      prev1=""
      prev2=""
      alignment1.reverse.zip(alignment2.reverse).each_with_index do |pair, index|
        k1, k2 = pair
        if k1=="-" then this1="-" else this1="*" end
        if k2=="-" then this2="-" else this2="*" end
        if this1!=prev1
          if k1=="-"
            a[:gaps]+=1
            prev1=this1
          else
            a[:exons]+=1
            prev1=this1
          end
        end
        if this2!=prev2
          if k2=="-"
            b[:gaps]+=1
            prev2=this2
          else
            b[:exons]+=1
            prev2=this2
          end
        end
        if index==0 and k1=="-"
          a[:start_gap]=1
        end
        if index==alignment1.length-1 and k1=="-"
          a[:end_gap]=1
        end

        if index==0 and k2=="-"
          b[:start_gap]=1
        end
        if index==alignment2.length-1 and k2=="-"
          b[:end_gap]=1
        end
        if k1==k2
          a[:overlaps]+=1
          b[:overlaps]+=1
        end

      end

      # p a
      # p b

      type = -1
      if a[:overlaps]==0
        type = 5
      elsif b[:start_gap]==1 and b[:end_gap]==1 and b[:gaps]==2 and b[:exons]==1 and
         a[:exons]==1 and a[:gaps]==0
        type = 1
      elsif a[:start_gap]==1 and a[:end_gap]==1 and a[:gaps]==2 and a[:exons]==1 and
         b[:exons]==1 and b[:gaps]==0
        type = 2
      elsif a[:exons]==1 and a[:gaps]==0 and b[:exons]==1 and b[:gaps]==1 and
            (b[:start_gap]==1 or b[:end_gap]==1)
        type = 1
      elsif b[:exons]==1 and b[:gaps]==0 and a[:exons]==1 and a[:gaps]==1 and
            (a[:start_gap]==1 or a[:end_gap]==1)
        type = 2
      elsif a[:exons]==1 and a[:gaps]==1 and a[:end_gap]==1 and
            b[:exons]==1 and b[:gaps]==1 and b[:start_gap]==1
        type = 3
      elsif a[:exons]==1 and a[:gaps]==1 and a[:start_gap]==1 and
            b[:exons]==1 and b[:gaps]==1 and b[:end_gap]==1
        type = 4
      else
        type = 5
      end

      consensus = ""
      if type < 0
        puts "something went really wrong"
      elsif type <= 4
        alignment1.reverse.zip(alignment2.reverse).each_with_index do |pair, index|
          if index==0
            if pair[0]=="-"
              consensus << pair[1]
            elsif pair[1]=="-"
              consensus << pair[0]
            else
              consensus << pair[0]
            end
          else
            if pair[0]=="-"
              consensus << pair[1][-1]
            elsif pair[1]=="-"
              consensus << pair[0][-1]
            else
              consensus << pair[0][-1]
            end
          end
        end
      end

      return [type, consensus]
    end

    def max3(a, b, c)
      if a>=b and a>=c
        return [a, "M"]
      elsif b>=a and b>=c
        return [b, "X"]
      elsif c>=a and c>=b
        return [c, "Y"]
      end
    end

    def initialize_matrix s1, s2
      matrix = []
      0.upto(s1) do |x|
        matrix[x] ||= []
        matrix[x][0] = { :score_m => -1e6,
                         :score_x => -@gap_open-(x-1)*@gap_extend,
                         :score_y => -1e6,
                         :arrow_m => "X",
                         :arrow_x => "X",
                         :arrow_y => "X"} # global
      end
      0.upto(s2) do |y|
        matrix[0] ||= []
        matrix[0][y] = { :score_m => -1e6,
                         :score_x => -1e6,
                         :score_y => -@gap_open-(y-1)*@gap_extend,
                         :arrow_m => "Y",
                         :arrow_x => "Y",
                         :arrow_y => "Y"} # global
      end
      matrix[0][0] = {:score_m => 0,
                      :score_x => 0,
                      :score_y => 0,
                      :arrow_m => "E",
                      :arrow_x => "E",
                      :arrow_y => "E"}
      return matrix
    end

    def kmerise id, seq
      if seq.nil?
        abort "sequence #{id} is nil!!!"
      end
      list = []
      (0..seq.length-@kmer).each do |i|
        kmer = (seq[i..(i+@kmer-1)]).upcase
        list << kmer
      end
      return list
    end

  end

end