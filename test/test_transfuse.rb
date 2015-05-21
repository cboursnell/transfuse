#!/usr/bin/env	ruby

require 'helper'
require 'tmpdir'

class TestTransfuse < Test::Unit::TestCase

  context 'transfuse' do

    setup do
      @fuser = Transfuse::Transfuse.new 4
    end

    teardown do
    end

    should 'check for existence of files' do
      list = []
      list << File.join(File.dirname(__FILE__), 'data', 'assembly1.fasta')
      list << File.join(File.dirname(__FILE__), 'data', 'assembly2.fasta')
      files = @fuser.check_files list.join(",")
      assert_equal 2, files.length, "length"
    end

    should "concatenate two files" do
      list = []
      list << File.join(File.dirname(__FILE__), 'data', 'assembly1.fasta')
      list << File.join(File.dirname(__FILE__), 'data', 'assembly2.fasta')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir(tmpdir) do
          output = @fuser.concatenate list
          assert File.exist?(output)
          lines = `wc -l #{output}`
          assert_equal 1000, lines.split(/\s+/)[0].to_i
        end
      end
    end

    # should "cluster fasta file" do
    #   Dir.mktmpdir do |tmpdir|
    #     Dir.chdir(tmpdir) do
    #       file = File.join(File.dirname(__FILE__), 'data', 'assembly1.fasta')
    #       hash = @fuser.cluster file
    #       assert_equal 250, hash.size, "output size"
    #     end
    #   end
    # end

    should "load scores from transrate output" do
      files = []
      files << File.join(File.dirname(__FILE__), 'data', 'contig_scores1.csv')
      hash = @fuser.load_scores files
      assert_equal 99, hash.size
    end

    should "select contigs" do
      clusters = {"0" => ["contig1", "contig2"], "1" => ["contig3", "contig4"]}
      scores = { "contig1" => 0.2,
                 "contig2" => 0.3,
                 "contig3" => 0.4,
                 "contig4" => 0.2 }
      best = @fuser.select_contigs clusters, scores
      assert_equal 2, best.size
      assert_equal "contig2", best[0]
      assert_equal "contig3", best[1]
    end

    should "output contigs" do
      best = ["soap_contig173359", "oases_contig80246"]
      file = File.join(File.dirname(__FILE__), 'data', 'assembly1.fasta')
      Dir.mktmpdir do |tmpdir|
        Dir.chdir(tmpdir) do
          @fuser.output_contigs best, file, "out"
        end
      end
    end

  end
end
