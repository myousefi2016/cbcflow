# For only selected schemes and problems:
#time py.test -s --testsize debug --schemes IPCS_Stable u_degree=2 --type update --problems Beltrami LidDrivenCavity
time py.test -s --testsize all --type update
