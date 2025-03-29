import numpy as np
import matplotlib.pyplot as plt
import random as rnd


def covid_grid(N):
    grid = np.zeros((N, N), dtype = int)
    return grid


def range_influence_prob(xi_r, R0):
    """
    probability for the range of influence for the day
    return as int to match indices
    """
    return -R0 * np.log(1-xi_r)


def mobility_prob(xi_m, T):
    """
    mobility probability. Returns m. Note, number
    of interactions given by m*sigma, where sigma is constant
    Returns
    """
    # still may be wrong i dont know
    return - T * np.log(1 - xi_m)


def select_random_candidates(n_infect, rzi):
    """
    n_infect : number of infected (number of random spots to be chosen)
    rzi      : zone of influence 

    returns an array with each entry being random i,j coordinates
    """
    indx_rand_list = np.zeros((n_infect, 2), dtype=int)

    for indx in range(n_infect):
        ri = int(np.rint(rnd.uniform(-rzi, rzi)))
        rj = int(np.rint(rnd.uniform(-rzi, rzi)))
        indx_rand_list[indx, 0] = ri
        indx_rand_list[indx, 1] = rj
    return indx_rand_list


def init_infect(grid, N_infect, infect_range):
    """
    Make sure infect range is half the length of the grid
    """
    candidate_list = select_random_candidates(N_infect, infect_range)
    for i in range(N_infect):
        ix, iy = candidate_list[i, 0], candidate_list[i, 1]
        grid[ix, iy] = 1.
    return grid


def grid_scanner(grid, time_grid, R0, T, sigma, pi, pd, bed_grid, p_icu, vuln_grid):
	"""
	Takes a grid and puts indices of all infecteds in an array. 
	Calculates range of influence and numbers infected around that cell.
	Repeats for all initially infected cells
	"""
	# Find infected cells
	inf_indices = np.where(grid == 1)
	inf_indices = np.vstack((inf_indices[0], inf_indices[1])).T # n by 2 grid

	# Initialize the counters for recovery or death
	recover_count = 0
	death_count = 0
	icu_count = 0 # Daily icu count init

	for infected in range(len(inf_indices)):
		# current infected cell in loop
		inf_x, inf_y = inf_indices[infected, 0], inf_indices[infected, 1]

		# Range of influence for this cell only
		Rzi = range_influence_prob(np.random.uniform(), R0)

		# Need the number of cells interacted with using mobi prob
		m = mobility_prob(np.random.uniform(), T)
		n = 1
		N_infect = int(np.rint(m * n * sigma)) 


		for candidates in range(N_infect):

			# new indices with respect to the zone of influence.
			# Rounded and converted to nearest integer
			new_inf_x = int(np.rint(rnd.uniform(-Rzi, Rzi)))	
			new_inf_y = int(np.rint(rnd.uniform(-Rzi, Rzi)))	

			# New candidate cell for infection. Need to run probability tho
			new_cell_x, new_cell_y = inf_x + new_inf_x, inf_y + new_inf_y

			# If statement is to stop index rollover.
			# Try/except/pass used to skip over out of range index errors.
			if new_cell_x >= 0 and new_cell_y >= 0:
				try:
					# Check if new cell is already infected. Then run probability
					if grid[new_cell_x, new_cell_y] == False:
						pi_check = np.random.uniform()
						if pi_check < pi:
							grid[new_cell_x, new_cell_y] = 1
							time_grid[new_cell_x, new_cell_y] = 21

						# Check to see if they're admitted to ICU
						if icu_check == True:
							if bed_grid[new_cell_x, new_cell_y] == False: # ie, not already in ICU
								xi_icu = np.random.uniform()
								if xi_icu <= p_icu: 
									bed_grid[new_cell_x, new_cell_y] = 1
									# Set to 2 so it cant reinfect
									grid[new_cell_x, new_cell_y] = 2 # TESTING REMOVING FROM GRID
									icu_count += 1

									# Should we take the cell out of the game?

				except:
					pass
		
	return grid, time_grid, bed_grid, recover_count, death_count, icu_count


def full_cycle(N_grid, R0, T, sigma, pi, pd, p_icu, pd_vuln, vuln_rate, days, lockdown = False, icu_check = False, vuln_check = False):
	# make grid and time grid
	grid = covid_grid(N_grid)
	time_grid = np.copy(grid)
	bed_grid = np.copy(grid)

	# Infect 1 space
	grid[N_grid//2, N_grid//2] = 1 #sets infected point to middle
	time_grid[N_grid//2, N_grid//2] = 21 

	# count lists for plotting
	death_count = []
	cum_death_count = []
	recover_count = []
	icu_count = []
	cum_icu_count = []
	day_array = []
	grid_history = [time_grid]

	# set total icu count
	total_icu_count = 0
	pd_original = pd # used for resetting pd after the vuln check.

	# Runs each day. Stores how long each cell has had covid in the time matrix
	for day in range(1, days+1):

		if lockdown:

			# Switches radius on day 21 for the lockdown
			if day == 21:
				R0 = 3

		grid, time_grid, bed_grid, recover_counter, death_counter, daily_icu_counter = grid_scanner(grid, 
																	                time_grid, 
																	                R0, 
																	                T, 
																	                sigma, 
																	                pi, 
																	                pd, 
																	                bed_grid, 
																	                p_icu, 
																	                icu_check)

		time_grid = time_grid - 1
	
		if 1 in time_grid:
			time_indices = np.where(time_grid == 1)
			time_indices = np.vstack((time_indices[0], time_indices[1])).T

			# Once the timer hits 1, find if they die and log it. reset grid points.
			for death_indx in range(len(time_indices)): 
				xi_death = np.random.uniform()

				# Frees up grid point.
				grid[time_indices[death_indx][0], time_indices[death_indx][1]] = 0
				bed_grid[time_indices[death_indx][0], time_indices[death_indx][1]] = 0

				# Counts death if they die, else they recover. Also checks ICU

				# Checks if the cell in question is vulnerable
				if vuln_check == True:
					xi_vuln = np.random.uniform()
					if xi_vuln <= vuln_rate:
						pd = pd_vuln

				if xi_death <= pd:
					death_counter += 1
					if bed_grid[time_indices[death_indx][0], time_indices[death_indx][1]] == True:
						daily_icu_counter -= 1
						total_icu_count -= 1

				else:
					recover_counter += 1
					# bed_grid[time_indices[death_indx][0], time_indices[death_indx][1]] = 0
				pd = pd_original # Resets the pd for death for the next cell

			total_icu_count = total_icu_count + daily_icu_counter
			cum_icu_count.append(total_icu_count)

			# I don't think we need recovery; i fucked up. will delete later.
			death_count.append(death_counter)
			recover_count.append(recover_counter)

		# All the counters and arrays for plotting deaths etc.	
		else:
			death_count.append(0)
			recover_count.append(0)
			cum_icu_count.append(total_icu_count)

		icu_count.append(daily_icu_counter)
		cum_death_count.append(np.sum(death_count))
		day_array.append(day)

		# Fixes the negative values in (unused) cells back to zero.
		time_grid[time_grid < 0] = 0
		grid_history.append(time_grid)
		
	return grid, grid_history, day_array, death_count, recover_count, cum_death_count, icu_count, cum_icu_count


if __name__ == '__main__':

	N = 300
	days = 200
	R0 = 8
	pi = 0.2
	pd = 0.3
	p_icu = 10**-2
	pd_vuln = 0.6
	vuln_rate = 0.1
	T = 20
	sigma = 0.1
	lockdown = False
	icu_check = True
	vuln_check = False


	covid, grid_history, days, death_count, recover_count, cum_death, icu_count, cum_icu = full_cycle(N, R0, T, sigma, pi, pd, p_icu, pd_vuln, vuln_rate, days, lockdown, icu_check, vuln_check)

	#plt.imshow(covid)#, cmap='Greys')
	#plt.colorbar()
	#plt.show()

	# days = np.arange(1, len(death_count) + 1)
	#np.savetxt(icu_count, cum_icu, "")
	# plt.plot(days, death_count)
	#plt.plot(days, cum_death)

	if icu_check == True:
		#Plot figure
		fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (10, 5))
		Nbeds = (N**2)/(1.97*1000)
		ax1.scatter(days, cum_icu, marker = 'd', color = "red", s = 0.1)
		ax1.plot(days, cum_icu, color = "red", lw = 0.1)
		ax1.hlines(y = Nbeds, xmin = 0, xmax = 200, linewidth = 1, color = 'black', ls = '--')
		
		ax2.bar(days, icu_count, color = "red", alpha = 0.6, width = 0.5, label = 'Daily ICU Rate') #marker='s',
		#ax2.plot(t[0:len(t)-(len(t)-len(icu_avg))], icu_avg, color = "black", label = '7 Day Average')
		ax2.hlines(y = Nbeds, xmin = 0, xmax = 200, linewidth = 1, color = 'black', ls = '--')

		# Set axes labels
		ax1.set_xlabel("Time (days)")
		ax2.set_xlabel("Time (days)")
		ax2.set_ylabel("Daily Hospitalization")
		ax1.set_ylabel("Cumulative Hospitalizations")
		#fig.suptitle("Simulated COVID-19 Hospitalizations")

		# Scale x axis limits to non-zero values
		#ax1.set_xlim(t[np.min(np.nonzero(icu_list))], t[np.max(np.nonzero(icu_list))])
		#ax1.set_xlim(t[np.min(np.nonzero(cum_icu))], t[np.max(np.nonzero(cum_icu))])
		#ax1.set_ylim()
		fig.savefig('icu_curve_plots2_1.pdf', dpi = 300)
		#plt.plot(days, cum_icu)
		#plt.plot(days, icu_count)
	
	plt.show()


