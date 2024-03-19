# Stage 1: Build STAR
FROM ubuntu:20.04

#Set environmental variables to use my tools on MGI
ENV WD=/usr/local/bin/ \
PATH=$PATH:/usr/local/bin/utilities:/gscuser/m.inkman/local/bin:/gscuser/m.inkman/downloads/subread_src/subread-1.6.1-source/bin:/gscuser/m.inkman/local/ncbi/sra-tools/bin:/gapp/noarch/bin:/gapp/x64linux/bin:/gapp/ia32linux/bin:/gsc/scripts/opt/genome/current/user/bin:/gsc/scripts/bin:/gsc/scripts/opt/lims/snapshots/current/bin:/gsc/bin:/gsc/java/bin:/opt/lsf9/9.1/linux2.6-glibc2.3-x86_64/bin:/opt/lsf9/9.1/linux2.6-glibc2.3-x86_64/etc:/gsc/scripts/gsc/maherlab:/gsc/pkg/oracle/10gR2/db_1/bin \
LD_LIBRARY_PATH=/usr/local/lib:/usr/lib:/gscuser/m.inkman/local/lib:/gscuser/m.inkman/local/ngs/lib64:/gscuser/m.inkman/local/ncbi/ncbi-vdb/lib64:/opt/lsf9/9.1/linux2.6-glibc2.3-x86_64/lib:/gapp/x64linux/lib:/gapp/ia32linux/lib:/gsc/lib:/gsc/pkg/oracle/10gR2/db_1/lib \
CLASSPATH=$CLASSPATH:/gscuser/m.inkman/local/ngs/jar/ngs-java.jar \
NGS_LIBDIR=/gscuser/m.inkman/local/ngs/lib64 \
DEBIAN_FRONTEND=noninteractive

WORKDIR $WD
ADD . $WD

RUN apt-get update --yes \
&& apt-get install wget build-essential gdb python-dev python3-pip python3-dev libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libkrb5-dev libnss3-dev ssh sudo less zip unzip git-all cmake libxml2-dev bedtools --yes \
&& apt-get clean all


# Install required python packages
RUN wget https://bootstrap.pypa.io/pip/2.7/get-pip.py \
&& python get-pip.py \
&& rm get-pip.py \
&& pip install numpy pysam matplotlib pandas scipy whichcraft \
&& pip3 install numpy pysam matplotlib pandas scipy whichcraft

# RUN sudo truncate -s 0 /etc/samba/smb.conf

# Install OpenSSL
ARG ver="openssl-1.1.1j"
RUN wget --no-check-certificate https://www.openssl.org/source/$ver.tar.gz \
&& tar xf $ver.tar.gz \
&& rm $ver.tar.gz \
&& cd $ver \
&& ./config --prefix=/usr --libdir=/lib/x86_64-linux-gnu \
&& make \
&& make install \
&& cd .. \
&& rm -rf $ver

# Install curl
RUN wget --no-check-certificate https://curl.haxx.se/download/curl-7.75.0.tar.gz \
&& tar xf curl-7.75.0.tar.gz \
&& rm curl-7.75.0.tar.gz \
&& cd curl-7.75.0 \
&& ./configure --prefix=/usr --libdir=/usr/lib/x86_64-linux-gnu --with-ssl --with-libssl-prefix=/usr/lib/x86_64-linux-gnu \
&& make \
&& make install \
&& cd .. \
&& rm -rf curl-7.75.0

# Install HDF5
RUN wget --no-check-certificate https://support.hdfgroup.org/ftp/HDF5/current18/src/hdf5-1.8.21.tar.gz \
&& tar xf hdf5-1.8.21.tar.gz \
&& rm hdf5-1.8.21.tar.gz \
&& cd hdf5-1.8.21 \
&& ./configure --prefix=/usr/local \
&& make \
&& make install \
&& cd .. \
&& rm -rf hdf5-1.8.21

ENV HDF5_DIR=/usr/local


# Set default values for ARGs
ARG ver_star=2.7.5b

# Install STAR
RUN wget --no-check-certificate https://github.com/alexdobin/STAR/archive/2.7.5b.tar.gz \
&& tar xf 2.7.5b.tar.gz \
&& mv STAR-2.7.5b/bin/Linux_x86_64/STAR* /usr/local/bin \
&& rm -rf STAR-2.7.5b 2.7.5b.tar.gz

# Install samtools
RUN wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
&& tar xf samtools-1.9.tar.bz2 \
&& rm samtools-1.9.tar.bz2 \
&& cd samtools-1.9 \
&& ./configure --prefix=/usr/local \
&& make all all-htslib \
&& make install install-htslib \
&& cd .. \
&& rm -rf samtools-1.9

# Set the working directory
WORKDIR $WD

# Copy only the necessary binaries from the previous stages
# COPY --from=star-builder /usr/local/bin/STAR* /usr/local/bin/
# COPY --from=samtools-builder /usr/local/bin/samtools /usr/local/bin/

# RUN apt-get update && apt-get install -y git

# Clone the repository
RUN git clone https://github.com/zacheliason/hla-em.git

# Set the working directory to the cloned repository
WORKDIR /usr/local/bin/hla-em

# Modify the Python script in-place
RUN sed -i -e '1i#!/usr/bin/env python3\n' HLA_EM.py

# Make the script executable
RUN chmod +x HLA_EM.py

# Define default command
CMD ["-h"]

# Set the entry point
ENTRYPOINT ["/usr/local/bin/hla-em/HLA_EM.py"]

