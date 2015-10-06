#!/usr/bin/env  ruby

require 'helper'
require 'tmpdir'

class TestCluster < Test::Unit::TestCase

  context 'cluster' do

    setup do
      @cluster = Transfuse::Cluster.new 4, true, 1.0
    end

    teardown do
    end

    should 'initialize' do
      assert @cluster
    end

    should 'generate vsearch command' do
      output = @cluster.generate_vsearch_command "assembly1.fasta", "output.txt", "output.msa"
      a = "vsearch --cluster_fast assembly1.fasta --id 1.0 "
      a << "--iddef 0 --qmask none --strand both --uc output.txt "
      a << "--msaout output.msa --threads 4"
      b = output.split(" ")
      b[0] = File.basename(b[0])
      output = b.join(" ")
      assert_equal a, output, "cmd"
    end

  end
end
