#!/usr/bin/env  ruby

require 'helper'
require 'tmpdir'

class TestCluster < Test::Unit::TestCase

  context 'cluster' do

    setup do
      @cluster = Transfuse::Cluster.new 4
    end

    teardown do
    end

    should 'initialize' do
      assert @cluster
    end

    should 'generate command' do
      output = @cluster.generate_command "assembly1.fasta", "output.txt"
      a = "vsearch --cluster_fast assembly1.fasta --id 1.00 "
      a << "--strand both --uc output.txt --threads 4"
      b = output.split(" ")
      b[0] = File.basename(b[0])
      output = b.join(" ")
      assert_equal a, output, "cmd"
    end

  end
end
