from matplotlib import pyplot as plt

plt.figure(1, figsize=(10,8), dpi=150)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.subplot(1,2,1)
plt.contourf(base.x_grid, base.y_grid, base_bamber)

plt.subplot(1,2,2)
plt.contourf(base.x_grid, base.y_grid, base.topg[:,:])

plt.show()

