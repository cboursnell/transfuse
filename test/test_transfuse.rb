#!/usr/bin/env	ruby

require 'helper'
require 'tmpdir'

class TestTransfuse < Test::Unit::TestCase

  context 'transfuse' do

    setup do
      @fuser = Transfuse::Transfuse.new 4, true
    end

    teardown do
    end

    should '1 check for existence of files' do
      list = []
      list << File.join(File.dirname(__FILE__), 'data', 'assembly1.fasta')
      list << File.join(File.dirname(__FILE__), 'data', 'assembly2.fasta')
      files = @fuser.check_files list.join(",")
      assert_equal 2, files.length, "length"
    end

    should "2 concatenate two files" do
      list = []
      list << File.join(File.dirname(__FILE__), 'data', 'assembly1.fasta')
      list << File.join(File.dirname(__FILE__), 'data', 'assembly2.fasta')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir(tmpdir) do
          output = @fuser.concatenate list
          assert File.exist?(output)
          lines = `wc -l #{output}`
          assert_equal 1000, lines.chomp.split.first.to_i
        end
      end
    end

    should "3 cluster fasta file" do
      # Dir.mktmpdir do |tmpdir|
      tmpdir = Dir.mktmpdir
        Dir.chdir(tmpdir) do
          file = File.join(File.dirname(__FILE__), 'data', 'assembly1.fasta')
          hash = @fuser.cluster file
          assert_equal 250, hash.size, "output size"
        end
      # end
    end

    should "4 load scores from transrate output" do
      files = []
      files << File.join(File.dirname(__FILE__), 'data', 'contig_scores1.csv')
      hash = @fuser.load_scores files
      assert_equal 99, hash.size
    end

    should "5 run transrate on assembly files with reads" do
      files = []
      left = []
      right = []
      files << File.join(File.dirname(__FILE__), 'data', 'assembly3.fasta')
      left << File.join(File.dirname(__FILE__), 'data', 'left.fq')
      right << File.join(File.dirname(__FILE__), 'data', 'right.fq')
      # Dir.mktmpdir do |tmpdir|
      tmpdir = Dir.mktmpdir
        Dir.chdir(tmpdir) do
          scores = @fuser.transrate files, left, right
          assert_equal 100, scores.size, "scores size"
        end
      # end
    end

    should "6 filter contigs" do
      files = []
      left = []
      right = []
      files << File.join(File.dirname(__FILE__), 'data', 'assembly1.fasta')
      left << File.join(File.dirname(__FILE__), 'data', 'left.fq')
      right << File.join(File.dirname(__FILE__), 'data', 'right.fq')
      # Dir.mktmpdir do |tmpdir|
      tmpdir = Dir.mktmpdir
        Dir.chdir(tmpdir) do
          scores = @fuser.transrate files, left, right
          scores.each do |contig, score|
            # puts "#{contig}\t#{score}"
          end
          new_list = @fuser.filter files, scores
          assert_equal 1, new_list.length
          cmd = "grep -c \">\" #{new_list.first}"
          assert_equal 1, `#{cmd}`.chomp.split.first.to_i, "number of contigs"
        end
      # end

    end

    should "7 get consensus of clusters" do

    end

    should "8 not fail when there are duplicated kmers in the input sequences" do

    end

  end
end
