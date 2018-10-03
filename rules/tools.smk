rule get_fasttree:
    input:
        HTTP.remote(
            'www.microbesonline.org/fasttree/FastTree',
            insecure=True, keep_local=True)
    output:
        'bin/FastTree'
    shell:
        'mv {input} {output} && chmod +x {output}'

rule get_mafft:
    input:
        HTTP.remote(
            'mafft.cbrc.jp/alignment/software/mafft-7.407-linux.tgz',
            insecure=True, allow_redirects=True)
    output:
        'bin/mafft-linux64/mafft.bat'
    shell:
        'tar xfvz {input} -C bin'

rule get_notung:
    input:
        HTTP.remote(
            'goby.compbio.cs.cmu.edu/Notung/Notung-2.9.zip',
            insecure=True)
    output:
        'bin/Notung-2.9/Notung-2.9.jar'
    shell:
        'unzip {input} -d bin'