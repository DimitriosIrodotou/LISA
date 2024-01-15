import re
import glob
import h5py
import numpy as np

# # # # # # # # # # Edit the variables below as desired # # # # # # # # # #
# Boolean variables to read in new log, subfind, or subhalo data (if False, load existing numpy files). #
read_log_data = False  # Read in black hole log files, extract, and save information during black hole mergers.
read_subfind_data = False  # Read in subfind data to extract galactic properties during black hole mergers.
avoid_overwriting = False  # If data for a region has been already saved, do not overwrite it and move to the next region.

# Path to save the data #
save_data_path = '/cosma7/data/dp004/dc-irod1/FLARES/LISA_data/'
# # # # # # # # # # Edit the variables above as desired # # # # # # # # # #

# Variables that contain info for all available FLARES data. Edit if a selection of redshifts and regions is needed. #
data_path = '/cosma7old/data/dp004/FLARES/FLARES-1/flares_'
particle_data_path = '/cosma7old/data/dp004/FLARES/FLARES-1/flares_'
redshifts = [4.77, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00]
tags = ['011_z004p770', '010_z005p000', '009_z006p000', '008_z007p000', '007_z008p000', '006_z009p000', '005_z010p000',
        '004_z011p000', '003_z012p000', '002_z013p000', '001_z014p000', '000_z015p000']
regions = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17',
           '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35',
           '36', '37', '38', '39']

# In the following block of code load, search, and extract from all black hole details files black hole properties
# written down at the time of black hole mergers (indicated by the word 'swallows' in the log files).
if read_log_data is True:
    # Loop over all regions. #
    for region in regions:
        print("Analysing blackhole_details for region", region)

        # Avoid overwriting existing black hole log file data. #
        if avoid_overwriting is True:
            # Check if a region's data already exists, if yes continue to the next region. #
            region_names = glob.glob(save_data_path + 'times_*')
            region_names = [re.split('_|.npy', name)[2] for name in region_names]
            if str(region) in region_names:
                continue

        # Path to the FLARES 'blackhole_details'. #
        local_data_path = data_path + region + '/data/blackhole_details/'

        # Declare arrays to store the information from the 'blackhole_details_' text files. #
        times, ids_primary, ids_secondary, masses_primary, masses_secondary = [], [], [], [], []

        # Find the files named 'blackhole_details_' and store their names. #
        blackhole_file_names = glob.glob(local_data_path + 'blackhole_details_*')
        blackhole_file_names = [re.split('/|./', name)[-1] for name in blackhole_file_names]

        # Loop over each 'blackhole_details_' text file. #
        for file_name in blackhole_file_names:
            with open(local_data_path + file_name) as file:
                # Loop over each line of each 'blackhole_details_' text file. #
                for line in file:
                    # Remove empty lines and split the line into words #
                    line = line.strip()  # Remove '\n' at end of line.
                    words = line.split(" ")

                    # Search for the word 'swallows' and get the relevant information from that line. #
                    if 'swallows' in words:
                        time = re.split('[=:]', words[1])[1]
                        id_primary = re.split('=|', words[2])[1]
                        id_secondary = words[4]
                        mass_primary = re.split(r'[(,]', words[5])[1]
                        mass_secondary = re.split(r'[),]', words[6])[0]

                        # Append the relevant information from each line for all files into arrays for each region. #
                        times.append(time)
                        ids_primary.append(id_primary)
                        ids_secondary.append(id_secondary)
                        masses_primary.append(mass_primary)
                        masses_secondary.append(mass_secondary)

            file.close()  # Close the opened file.

        # Save data in numpy arrays #
        np.save(save_data_path + 'times_' + str(region), times)
        np.save(save_data_path + 'ids_primary_' + str(region), ids_primary)
        np.save(save_data_path + 'ids_secondary_' + str(region), ids_secondary)
        np.save(save_data_path + 'masses_primary_' + str(region), masses_primary)
        np.save(save_data_path + 'masses_secondary_' + str(region), masses_secondary)

# In the following block of code load the black hole properties saved above and match the times of black hole mergers
# to snapshot files. Then extract the required properties (e.g. group and subgroup numbers) from the snapshot to be
# able to match black holes to their host galaxies.
if read_subfind_data is True:
    # Loop over all regions. #
    for region in regions:
        # Avoid analysing non-existing regions. #
        try:
            times = np.load(save_data_path + 'times_' + str(region) + '.npy')
        except FileNotFoundError:
            continue
        if len(times) == 0: continue

        # Load the remaining data saved in the read-in block of code. #
        ids_primary = np.load(save_data_path + 'ids_primary_' + str(region) + '.npy')
        ids_secondary = np.load(save_data_path + 'ids_secondary_' + str(region) + '.npy')
        masses_primary = np.load(save_data_path + 'masses_primary_' + str(region) + '.npy')
        masses_secondary = np.load(save_data_path + 'masses_secondary_' + str(region) + '.npy')

        # Convert the data into arrays of floats. #
        times = np.array(times, dtype=float)
        ids_primary = np.array(ids_primary, dtype=int)
        ids_secondary = np.array(ids_secondary, dtype=int)
        masses_primary = np.array(masses_primary, dtype=float)
        masses_secondary = np.array(masses_secondary, dtype=float)

        # Mask the data to prevent cases where the secondary black hole has been turned into a ghost particles. #
        mask_valid, = np.where((masses_primary > 0) & (masses_secondary > 0))
        times = times[mask_valid]
        ids_primary = ids_primary[mask_valid]
        ids_secondary = ids_secondary[mask_valid]
        masses_primary = masses_primary[mask_valid]
        masses_secondary = masses_secondary[mask_valid]

        # Round redshifts up/down to the closest integer, but exclude values that need to be assigned to redshift 4.77 #
        bh_redshifts = 1 / times - 1  # Convert times (i.e. scale factors) into redshift.
        mask_low_z, = np.where(bh_redshifts <= 4.885)  # 4.885 = 4.77 + (5.00 - 4.77) / 2
        stellar_mass_redshifts = np.round(bh_redshifts, decimals=0)
        stellar_mass_redshifts[mask_low_z] = 4.77

        # Declare an array to store the redshifts converted into tags. #
        snapshots = np.array(["000_z000p000" for x in range(len(stellar_mass_redshifts))])

        # Use the rounded up/down redshifts to load snapshot data and extract particle properties. #
        for i, (redshift, tag) in enumerate(zip(redshifts, tags)):
            mask_redshift, = np.where(redshift == stellar_mass_redshifts)
            snapshots[mask_redshift] = tags[i]

        for redshift in redshifts:
            # Avoid overwriting existing subfind data. #
            if avoid_overwriting is True:
                # Check if a region's data already exists, if yes continue to the next region. #
                redshift_names = glob.glob(save_data_path + 'ids_' + str(region) + '_' + str(redshift) + '.npy')
                redshift_names = [re.split('_|.npy', name)[3] for name in redshift_names]
                if str(redshift) in redshift_names:
                    continue

            # Declare arrays to store the information from the 'eagle_subfind_particles_' text files. #
            ids, masses, group_numbers, subgroup_numbers, particledata_files = [], [], [], [], []

            # Avoid analysing non-existing redshifts. #
            if redshift not in stellar_mass_redshifts: continue

            mask_redshift, = np.where(redshift == stellar_mass_redshifts)
            print("Found", len(mask_redshift), "black hole mergers at redshift", redshift, "in region", region)

            # Path to the 'particledata'. #
            local_particle_data_path = particle_data_path + region + '/data/particledata_' \
                                       + snapshots[mask_redshift][0] + '/eagle_subfind_particles_' \
                                       + snapshots[mask_redshift][0] + '.'

            # Find the files named 'blackhole_details_' and store their names. #
            particledata_file_names = glob.glob(local_particle_data_path + '*.hdf5')
            particledata_file_names = [re.split(r'\.|\.hdf5', name)[-2] for name in particledata_file_names]

            # Loop over each 'eagle_subfind_particles_' file. #
            for file_name in particledata_file_names:
                # Open the file and get the black hole IDs. #
                file = h5py.File(local_particle_data_path + file_name + '.hdf5', 'r')
                # Avoid loading files without black holes. #
                try:
                    pt5_ids = file['PartType5/ParticleIDs']
                except KeyError:
                    continue
                pt5_masses = np.array(file['PartType5/Mass'])
                pt5_group_number = np.array(file['PartType5/GroupNumber'])
                pt5_subgroup_number = np.array(file['PartType5/SubGroupNumber'])

                # Append the relevant information from each file for each redshift into arrays for each region. #
                for entries in range(len(pt5_group_number)):
                    particledata_files.append(file_name)
                ids.append(pt5_ids)
                masses.append(pt5_masses)
                group_numbers.append(pt5_group_number)
                subgroup_numbers.append(pt5_subgroup_number)

            # Save data in numpy arrays #
            np.save(save_data_path + 'ids_' + str(region) + '_' + str(redshift), np.concatenate(ids))
            np.save(save_data_path + 'files_' + str(region) + '_' + str(redshift), particledata_files)
            np.save(save_data_path + 'masses_' + str(region) + '_' + str(redshift), np.concatenate(masses))
            np.save(save_data_path + 'group_numbers_' + str(region) + '_' + str(redshift),
                    np.concatenate(group_numbers))
            np.save(save_data_path + 'subgroup_numbers_' + str(region) + '_' + str(redshift),
                    np.concatenate(subgroup_numbers))