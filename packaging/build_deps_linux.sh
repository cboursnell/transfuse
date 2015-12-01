echo "don't use this"

export PATH=~/gcc-4.8.4/bin:~/cmake-3.1.3-Linux-x86_64/bin:$PATH

cd ~

# vsearch
wget https://github.com/torognes/vsearch/releases/download/v1.8.0/vsearch-1.8.0-linux-x86_64.tar.gz
tar xzf vsearch-1.8.0-linux-x86_64.tar.gz
cd vsearch-1.8.0-linux-x86_64
cd bin
tar zcvf vsearch.tar.gz vsearch
## package it up
cp vsearch.tar.gz /vagrant/



