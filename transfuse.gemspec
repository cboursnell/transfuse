
require File.expand_path('../lib/transfuse/version', __FILE__)

Gem::Specification.new do |gem|
  gem.name        = 'transfuse'
  gem.version     = Transfuse::VERSION::STRING.dup
  gem.date        = '2015-05-22'
  gem.summary     = "Merge assemblies"
  gem.description = "See summary"
  gem.authors     = ["Richard Smith-Unna", "Chris Boursnell"]
  gem.email       = ['rds45@cam.ac.uk', 'cmb211@cam.ac.uk']
  gem.files       = `git ls-files`.split("\n")
  gem.executables = ["transfuse"]
  gem.require_paths = %w( lib )
  gem.homepage    = 'https://github.com/cboursnell/transfuse'
  gem.license     = 'MIT'

  gem.add_dependency 'trollop', '~> 2.0'
  gem.add_dependency 'bio', '~> 1.4', '>= 1.4.3'
  gem.add_dependency 'fixwhich', '~> 1.0', '>= 1.0.2'
  gem.add_dependency 'bindeps', '~> 1.0', '>= 1.0.1'
  gem.add_dependency 'transrate', '~> 1.0', '>= 1.0.1'

  gem.add_development_dependency 'rake', '~> 10.3', '>= 10.3.2'
  gem.add_development_dependency 'turn', '~> 0.9', '>= 0.9.7'
  gem.add_development_dependency 'simplecov', '~> 0.8', '>= 0.8.2'
  gem.add_development_dependency 'shoulda-context', '~> 1.2', '>= 1.2.1'
  gem.add_development_dependency 'coveralls', '~> 0.7'
end
