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
  task :linux => [:bundle_install, "packaging/packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-linux-x86_64.tar.gz"] do
    create_package("linux-x86_64")
  end

  desc "Package your app for OS X"
  task :osx => [:bundle_install, "packaging/packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-osx.tar.gz"] do
    create_package("osx")
  end

  file "packaging/packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-linux-x86_64.tar.gz" do
    download_runtime("linux-x86_64")
  end

  file "packaging/packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-osx.tar.gz" do
    download_runtime("osx")
  end

  desc "Install gems to local directory"
  task :bundle_install do
    if RUBY_VERSION !~ /^2\.2\./
      abort "You can only 'bundle install' using Ruby 2.2, because that's what Traveling Ruby uses."
    end
    Bundler.with_clean_env do
      sh "env BUNDLE_IGNORE_CONFIG=1 bundle install --path packaging/vendor"
    end
    sh "rm -f packaging/vendor/*/*/cache/*"
  end
end

def create_package(target)
  package_pref = "#{PACKAGE_NAME}-#{VERSION}-#{target}"
  package_dir = File.join("packaging", package_pref)
  sh "rm -rf #{package_dir}"
  sh "mkdir -p #{package_dir}/lib/app"

  # copy things from gem to package
  sh "cp -r bin #{package_dir}/lib/" # bin
  sh "cp -r lib #{package_dir}/lib/" # lib
  sh "cp -r deps #{package_dir}/lib/" # deps
  sh "cp files.txt #{package_dir}/lib/" # deps
  sh "cp Gemfile* #{package_dir}/lib/" # Gemfiles
  sh "cp *gemspec #{package_dir}/lib/" # gemspec

  # download travelling ruby
  sh "mkdir #{package_dir}/lib/ruby"
  sh "tar -xzf packaging/packaging/traveling-ruby-#{TRAVELING_RUBY_VERSION}-#{target}.tar.gz -C #{package_dir}/lib/ruby"
  sh "cp packaging/transfuse #{package_dir}/transfuse"
  sh "cp -pR packaging/vendor/* #{package_dir}/lib/"
  sh "cp Gemfile Gemfile.lock #{package_dir}/lib/"

  # install binary dependencies
  sh "mkdir -p #{package_dir}/bin"
  sh "mkdir -p #{package_dir}/lib"
  sh "wget -nc https://github.com/Blahah/snap/releases/download/v1.0beta.18/snap_v1.0beta.18_linux.tar.gz"
  sh "wget -nc https://github.com/Blahah/transrate-tools/releases/download/v1.0.0/bam-read_v1.0.0_linux.tar.gz"
  sh "wget -nc https://github.com/COMBINE-lab/salmon/releases/download/v0.4.2/SalmonBeta-0.4.2_DebianSqueeze.tar.gz"
  sh "wget -nc https://github.com/torognes/vsearch/releases/download/v1.9.5/vsearch-1.9.5-linux-x86_64.tar.gz"
  sh "find . -maxdepth 1 -name '*.tar.gz' -exec tar xzf '{}' \\;" # unpack
  sh "cp snap-aligner #{package_dir}/bin/."
  sh "cp bam-read #{package_dir}/bin/."
  sh "cp SalmonBeta-0.4.2_DebianSqueeze/bin/salmon #{package_dir}/bin/."
  sh "cp -r SalmonBeta-0.4.2_DebianSqueeze/lib/* #{package_dir}/lib/."
  sh "cp vsearch-1.9.5-linux-x86_64/bin/vsearch #{package_dir}/bin/."

  sh "cp packaging/libruby.* #{package_dir}/lib/."

  sh "mkdir #{package_dir}/lib/.bundle"
  sh "cp packaging/bundler-config #{package_dir}/lib/.bundle/config"
  # create package
  if !ENV['DIR_ONLY']
    sh "cd packaging && tar -czf #{package_pref}.tar.gz #{package_pref}"
    sh "rm -rf #{package_dir}"
  end
  sh "rm -rf packaging/vendor packaging/bindeps .bundle"

end

def download_runtime(target)
  sh "mkdir -p packaging/packaging &&" +
     "cd packaging/packaging && curl -L -O --fail " +
     "http://d6r77u77i8pq3.cloudfront.net/releases/traveling-ruby-#{TRAVELING_RUBY_VERSION}-#{target}.tar.gz"
end
