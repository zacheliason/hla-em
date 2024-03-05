# Stage 1: Build STAR
FROM ubuntu:20.04 AS star-builder

# Set default values for ARGs
ARG ver_star=2.7.5b

# Download and install STAR
RUN /bin/sh -c "apt-get update && apt-get install -y wget tar \
    && wget --no-check-certificate https://github.com/alexdobin/STAR/archive/${ver_star}.tar.gz \
    && tar xf ${ver_star}.tar.gz \
    && mv STAR-${ver_star}/bin/Linux_x86_64/STAR* /usr/local/bin \
    && rm -rf STAR-${ver_star} ${ver_star}.tar.gz"

RUN apt-get update && apt-get install -y \
    libc6 \
    libc6-dev \
    && rm -rf /var/lib/apt/lists/*

# Stage 2: Build samtools
FROM ubuntu:20.04 AS samtools-builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    tar \
    bzip2 \
    libncurses5-dev \ 
    zlib1g-dev \ 
    libbz2-dev \ 
    liblzma-dev \ 
    && rm -rf /var/lib/apt/lists/*

# Download samtools source
RUN wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2

# Extract samtools source
RUN tar xf samtools-1.9.tar.bz2

# Remove downloaded tarball
RUN rm samtools-1.9.tar.bz2

# Move to samtools directory
WORKDIR /samtools-1.9

# Configure samtools
RUN ./configure --prefix=/usr/local

# Build samtools
RUN make all all-htslib

# Install samtools
RUN make install install-htslib

# Debugging steps (comment out as needed)
RUN ls -la /usr/local/bin     # Check contents of /usr/local/bin
RUN which samtools            # Check the path of the installed samtools


# Stage 3: Install pip3 and whichcraft
FROM ubuntu:20.04

# Install pip3 and other dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

ENV PATH="/usr/bin/python3:${PATH}"

# Install whichcraft using pip3
RUN pip3 install whichcraft matplotlib
RUN pip3 install matplotlib

# Copy only the necessary binaries from the previous stages
COPY --from=star-builder /usr/local/bin/STAR* /usr/local/bin/
COPY --from=samtools-builder /usr/local/bin/samtools /usr/local/bin/

# Set the working directory
WORKDIR /usr/local/bin

# Entry point
CMD ["bash"]
