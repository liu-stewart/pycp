"""This module contains the script for the calculator."""
from __future__ import annotations
from pycp.method.params import setting
import pathlib


class Script():
    """This class used to store the run script."""

    def __init__(self,
                 queue_name: str = "mem128",
                 number_node: int = 1,
                 cpu_per_node: int = 0,
                 source_list=[],
                 module_list=["vasp/5.4.4"],
                 run_list=["mpirun -np $NP vasp_std > cal.out"],
                 custom_flags=[]
                 ):
        """Init this class."""
        self.queue_name = queue_name
        self.number_node = number_node
        if cpu_per_node == 0:
            cpu_per_node = setting["cpu_per_node"][self.queue_name]
        self.cpu_per_node = cpu_per_node
        self.source_list = source_list
        self.module_list = module_list
        self.run_list = run_list
        self.custom_flags = custom_flags

    def __str__(self) -> str:
        """Return the string representation of the script."""
        string = "#!/bin/sh\n"
        string += f"#BSUB -q {self.queue_name}\n"
        string += f"#BSUB -n {self.number_node*self.cpu_per_node}\n"
        string += f"#BSUB -e %J.err\n"
        string += f"#BSUB -o %J.out\n"
        string += f"#BSUB -R \"span[ptile={self.cpu_per_node}]\"\n"
        for flag in self.custom_flags:
            string = string + flag + "\n"
        string += "hostfile=`echo $LSB_DJOB_HOSTFILE`\n"
        string += "NP=`cat $hostfile | wc -l`\n"
        string += "cd $LS_SUBCWD\n"
        string += "\n"
        for source in self.source_list:
            string += f"source {source}\n"
        string += "\n"
        for module in self.module_list:
            string += f"module load {module}\n"
        string += "\n"
        for run in self.run_list:
            string += run + "\n"
        return string

    def write(self, to_file: pathlib.Path | str = "script.lsf"):
        """Write the script to a file."""
        if isinstance(to_file, str):
            to_file = pathlib.Path(to_file)
        if not to_file.parent.exists():
            to_file.parent.mkdir(parents=True)
        with open(to_file, "w") as file:
            file.write(str(self))
