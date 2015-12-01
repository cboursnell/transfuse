su vagrant

# Setup
sudo apt-get update -y
sudo apt-get install -y curl build-essential git zlib1g zlib1g-dev cmake

# Install a multi-user RVM install with the latest ruby
wget https://keybase.io/mpapis/key.asc
gpg --import key.asc

\curl -sSL https://get.rvm.io | bash -s stable --ruby
#source /usr/local/rvm/scripts/rvm
source /home/vagrant/.rvm/scripts/rvm

# Add user to RVM group
#sudo usermod -a -G rvm vagrant

# Install dependencies
gem install bundler rake


# Install transrate
git clone https://github.com/Blahah/transrate.git
cd transrate
gem build *spec
gem install *gem
rm *gem
cd ..

# Install transfuse
git clone https://github.com/cboursnell/transfuse.git
cd transfuse
sudo chmod -R 775 .


bundle exec rake package:linux
