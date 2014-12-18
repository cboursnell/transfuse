require 'rake/testtask'

Rake::TestTask.new do |t|
  t.libs << 'test'
end

Rake::TestTask.new do |t|
  t.name = :corset
  t.libs << 'test'
  t.test_files = ['test/test_corset.rb']
end

Rake::TestTask.new do |t|
  t.name = :cluster
  t.libs << 'test'
  t.test_files = ['test/test_cluster.rb']
end

desc "Run tests"
task :default => :test
