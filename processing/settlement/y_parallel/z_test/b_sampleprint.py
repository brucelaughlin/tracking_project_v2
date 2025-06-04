sample_start = 1000
sample_end = 1010

print("polygons_settled_per_particle")
print(polygons_settled_per_particle[sample_start:sample_end,:])
print("particle_settle_times_per_pld")
print(particle_settle_times_per_pld[sample_start:sample_end,:])
print("stagnant_settlers_per_particle_per_pld")
print(stagnant_settlers_per_particle_per_pld[sample_start:sample_end,:])
print("particle_distances_per_pld_settlers")
print(particle_distances_per_pld_settlers[sample_start:sample_end,:].astype(int))
print("particle_distances_per_pld_allReleases")
print(particle_distances_per_pld_allReleases[sample_start:sample_end,:].astype(int))
print("particle_distances_per_pld_settlers_noStagnation")
print(particle_distances_per_pld_settlers_noStagnation[sample_start:sample_end,:].astype(int))
