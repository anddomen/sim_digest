# sim_digest
Simulated restriction enzyme digest

# This software is Copyright (c) 2022 Oregon State University.
# All Rights Reserved. You may use, modify and redistribute
# this software under the terms of the GNU Affero General
# Public License version 3 published by the Free Software
# Foundation. If you need any other terms of use, please
# contact the author or advantage@oregonstate.edu

This script uses packages BiocManager and SimRAD to simulate sequences and then digest them with restriction enzymes. You can choose a single or double digest. There are currently 10 supported enzymes: ApeKI, EcoRI, MseI, TaqI, SbfI, PstI, AciI, AgeI, NcoI, NotI. If the enzyme you need isn't on this list, you can input the cutsites yourself. It will give you the sequenceable bases, genome coverage, and depth. It will also return a graph of the size selected fragments. Please note, if you have a mac with Apple Silicone, this won't work on your computer at the time of this publishing. 
