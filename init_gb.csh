#!/bin/csh

# --- Configuration Variables ---
set TARGET_DIR = bfe-gb
set SYS_NAME = "my_protein_ligand" # <-- NEW: This name is now configurable!
# -------------------------------

# 1. Create the directory
mkdir -p $TARGET_DIR

# 2. Copy the required input files (Assuming index.ndx is present for the gmx_MMPBSA command)
cp step3_input.pdb $TARGET_DIR/
cp step5_1.tpr $TARGET_DIR/
cp nopbc.xtc $TARGET_DIR/
cp -r toppar $TARGET_DIR/
cp topol.top $TARGET_DIR/
cp index.ndx $TARGET_DIR/ 

# 3. Create the mmgbsa.in input file inside the target directory
#    Note: The $SYS_NAME variable is substituted here by the C Shell.
cat << EOF > $TARGET_DIR/mmgbsa.in
&general
sys_name="$SYS_NAME",
startframe=900,
endframe=1000,
interval = 1,
temperature=300,

/

&gb
igb=5, saltcon=0.150, intdiel=4.0
/
EOF

# 4. Create the run.csh execution script inside the target directory
cat << EOF_RUN > $TARGET_DIR/run.csh
#!/bin/csh
# This script executes the gmx_MMPBSA calculation.

mpirun -np 16 gmx_MMPBSA -O -i mmgbsa.in -cs step5_1.tpr -ci index.ndx -cg 1 14 -ct nopbc.xtc -cp topol.top
EOF_RUN

# 5. Make the execution script runnable
chmod +x $TARGET_DIR/run.csh

echo "Setup complete. The directory '$TARGET_DIR' is ready with sys_name set to '$SYS_NAME'."
