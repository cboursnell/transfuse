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
