Building a new kraken1 / krakenuniq database
============================================

Prerequisites:

* gcc (e.g., gcc-9 or equivalent supporting C++17, `click here <https://github.com/seqan/seqan3#dependencies>`_ for details)
* dustmasker
* GNU parallel
* krakenuniq (krakenuniq-download has to be in PATH, even if you only want to use kraken1)
* `KMC3 <https://github.com/refresh-bio/KMC/releases/download/v3.1.1/KMC3.1.1.linux.tar.gz>`_
* perl (for the krakenuniq-download script)

Build:

::

    $ git clone --recursive https://github.com/cpockrandt/kraken-decontaminated-database.git
    $ mkdir kraken-decontaminated-database/build && cd kraken-decontaminated-database/build
    $ cmake .. -DCMAKE_BUILD_TYPE=Release
    $ make -j4

Before running the script, please check the path to KMC and the number of threads in the first lines of ``build.sh``.

::

    $ ./build.sh
