rule Relatedness_check:
    threads: config["threads"]
    resources: mem_mb=config["mem_mb"]
    shell: "echo running with {threads} threads"

