module Transfuse

  require 'csv'

  class Transrate

    def initialize threads
      @threads = threads
      gem_dir = Gem.loaded_specs['transrate'].full_gem_path
      @transrate = File.join(gem_dir, "bin", "transrate")
    end

    def run files, left, right
      score_files = []
      dir = "transrate"
      FileUtils.mkdir_p(dir)
      Dir.chdir(dir) do
        files.each do |fasta|
          name = File.basename(fasta, File.extname(fasta))
          FileUtils.mkdir_p(name)
          Dir.chdir(name) do
            cmd = "#{@transrate} "
            cmd << " --assembly #{fasta}"
            cmd << " --left #{left.join(",")}"
            cmd << " --right #{right.join(",")}"
            cmd << " --outfile transrate"
            cmd << " --threads #{@threads}"
            outfile = "transrate_#{File.basename(fasta)}_contigs.csv"

            rater = Cmd.new(cmd)
            rater.run outfile
            if rater.status.success?
              File.open("#{File.basename(fasta)}.log","wb") do |out|
                out.write rater.stdout
              end
            else
              abort "Something went wrong running transrate\n#{rater.stderr}"
            end
            score_files << File.expand_path(outfile)
          end
        end
      end
      return score_files

    end

  end

end