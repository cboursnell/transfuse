require 'rake/testtask'

Rake::TestTask.new do |t|
  t.libs << 'test'
end

Rake::TestTask.new do |t|
  t.name = :cluster
  t.libs << 'test'
  t.test_files = ['test/test_cluster.rb']
end

Rake::TestTask.new do |t|
  t.name = :fuse
  t.libs << 'test'
  t.test_files = ['test/test_transfuse.rb']
end

Rake::TestTask.new do |t|
  t.name = :cons
  t.libs << 'test'
  t.test_files = ['test/test_consensus.rb']
end

desc "Run tests"
task :default => :test

# packaging

require 'bundler/setup'

PACKAGE_NAME = "transfuse"
VERSION = "0.5.0"
TRAVELING_RUBY_VERSION = "20150210-2.2.0"

desc "Package your app"
task :package => ['package:linux', 'package:osx']

namespace :package do

  desc "Package your app for Linux x86_64"
  task :linux => [:bundle_install, "packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-linux-x86_64.tar.gz"] do
    create_package("linux-x86_64")
  end

  desc "Package your app for OS X"
  task :osx => [:bundle_install, "packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-osx.tar.gz"] do
    create_package("osx")
  end

  desc "Install gems to local directory"
  task :bundle_install do
    if RUBY_VERSION !~ /^2\.2\./
      abort "You can only 'bundle install' using Ruby 2.2, because that's what Traveling Ruby uses."
    end
    Bundler.with_clean_env do
      sh "env BUNDLE_IGNORE_CONFIG=1 bundle install --path packaging/vendor --without development"
    end
    sh "rm -f packaging/vendor/*/*/cache/*"
  end
end

file "packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-linux-x86_64.tar.gz" do
  download_runtime("linux-x86_64")
end

file "packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-osx.tar.gz" do
  download_runtime("osx")
end

def create_package(target)
  package_dir = "#{PACKAGE_NAME}-#{VERSION}-#{target}" #eg transfuse-0.5-linux
  sh "rm -rf #{package_dir}"
  sh "mkdir -p #{package_dir}/lib/app"

  # copy things from gem to package
  sh "cp bin  #{package_dir}/lib/app/" # bin
  sh "cp lib  #{package_dir}/lib/app/" # lib
  sh "cp deps  #{package_dir}/lib/app/" # deps
  sh "cp files.txt  #{package_dir}/lib/app/" # deps
  sh "cp Gemfile*  #{package_dir}/lib/app/" # Gemfiles
  sh "cp *gemspec  #{package_dir}/lib/app/" # gemspec
  #

  sh "mkdir #{package_dir}/lib/ruby"
  sh "tar -xzf packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-#{target}.tar.gz -C #{package_dir}/lib/ruby"
  sh "cp packaging/wrapper.sh #{package_dir}/#{PACKAGE_NAME}"
  sh "cp -pR packaging/vendor #{package_dir}/lib/"
  sh "cp Gemfile Gemfile.lock #{package_dir}/lib/vendor/"

  # install binary dependencies
  sh "mkdir -p packaging/downloads"
  sh "mkdir -p packaging/bindeps/#{target}"
  sd "cd packaging/downloads"
  sh "wget https://github.com/Blahah/snap/releases/download/v1.0beta.18/snap_v1.0beta.18_linux.tar.gz"
  sh "wget https://github.com/Blahah/transrate-tools/releases/download/v1.0.0/bam-read_v1.0.0_linux.tar.gz"
  sh "wget https://github.com/COMBINE-lab/salmon/releases/download/v0.4.2/SalmonBeta-0.4.2_DebianSqueeze.tar.gz"
  sh "wget https://github.com/torognes/vsearch/releases/download/v1.8.1/vsearch-1.8.1-linux-x86_64.tar.gz"
  sh "find . -maxdepth 1 -name '*.tar.gz' -exec tar xzf '{}' \\;" # unpack
  sh "cp snap-aligner ../bindeps/#{target}/."
  sh "cp bam-read ../bindeps/#{target}/."
  sh "cp -r SalmonBeta-0.4.2_DebianSqueeze/bin ../bindeps/#{target}/."
  sh "cp -r SalmonBeta-0.4.2_DebianSqueeze/lib ../bindeps/#{target}/."
  sh "cp vsearch-1.8.1-linux-x86_64/bin/vsearch ../bindeps/#{target}/."
  sh "cp -r packaging/bindeps/#{target}/{bin,lib} #{package_dir}/"

  sh "mkdir #{package_dir}/lib/vendor/.bundle"
  sh "cp packaging/bundler-config #{package_dir}/lib/vendor/.bundle/config"
  # create package
  if !ENV['DIR_ONLY']
    sh "tar -czf #{package_dir}.tar.gz #{package_dir}"
    sh "rm -rf #{package_dir}"
  end
end

def download_runtime(target)
  sh "cd packaging && curl -L -O --fail " +
    "http://d6r77u77i8pq3.cloudfront.net/releases/traveling-ruby-#{TRAVELING_RUBY_VERSION}-#{target}.tar.gz"
end

