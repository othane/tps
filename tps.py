"""
simluation of any number of teathered particles
"""

from pylab import *
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class particle(object):

	ke = 8.9875517873681764e9

	def __init__(self, q=None,m=None,x=None,y=None,z=None):
		self.q = self._rand_on_none_(q, xmin=10e-3, xmax=10.0) # charge in coulomb's
		self.m = self._rand_on_none_(m, xmin=0.50, xmax=1.0) # mass in kg
		self.pos = array([self._rand_on_none_(s, xmin=-1.0, xmax=1.0) for s in [x,y,z]])
		print "initial positions are ", self.pos, " with a norm of ", norm(self.pos)
		self.v = zeros(3) # initial velocity of the particles
		self._Fnew_ = zeros(3) # force on the particle
		self.F = self._Fnew_

	@staticmethod
	def _rand_on_none_(x, xmin=-1, xmax=1):
		w = xmax - xmin
		return w*rand(1)[0] + xmin if x == None else x

	def add(self, a):
		"""
		add the force from another particle (charge only)
		"""
		s = self.pos - a.pos
		sh = s / norm(s)
		self._Fnew_ += self.ke * (self.q * a.q) * sh / norm(s)**2

	@property
	def force(self):
		"""
		return the current force in Newtons acting on this particle
		"""
		return norm(self.F)

	@property
	def Fdir(self):
		fdir = self.pos - 0.6 * self.F / norm(self.F)
		return fdir

	def move(self, dt):
		"""
		Move the particle based on the old velocity and force from the other
		charges
		"""
		a = self._Fnew_ / self.m
		self.F = self._Fnew_
		self._Fnew_ = zeros(3) # set the force back to 0 for next update round
		self.pos += 0.5*a*dt**2 + self.v*dt
		self.v += a*dt

		# simulate the teather
		if norm(self.pos) > 1:
			self.pos *= 1.0 / norm(self.pos)

# stollen from http://stackoverflow.com/questions/22867620/putting-arrowheads-on-vectors-in-matplotlibs-3d-plot
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def set_3d_position(self, xs, ys, zs):
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

class tpsim(object):

	def __init__(self, N, ts = 0.1, to = 60.0, save=False):
		"""
		simulate N particles connected by a teather at a central point
		and update the simulation every ts seconds for up to to seconds
		"""
		self.ts = ts
		self.to = to
		self.particles = [particle(q=10e-6,m=1e-19) for n in range(N)]
		self.fig = figure()
		self.ax = self.fig.add_axes([0, 0, 1, 1], projection='3d')
		#self.ax.axis('off')
		self.colors = cm.Dark2(linspace(0, 1, len(self.particles)))
		self.pts = [self.ax.plot([], [], [], 'o', c=c)[0] for c in self.colors]
		self.lines = [self.ax.plot([], [], [], '-', c=c)[0] for c in self.colors]
		self.force_labels = [self.ax.text3D(p.pos[0]/2, p.pos[1]/2, p.pos[2]/2, "", color=c, size=10) for p,c in zip(self.particles, self.colors)]
		self.force_vectors = [self.ax.add_artist(Arrow3D([p.pos[0], p.Fdir[0]], [p.pos[1], p.Fdir[1]], [p.pos[2], p.Fdir[2]], axes=self.ax, mutation_scale=20, lw=0.5, arrowstyle="<-", color=c)) for p,c in zip(self.particles, self.colors)]
		self.ax.set_xlim(-2, 2)
		self.ax.set_ylim(-2, 2)
		self.ax.set_zlim(-2, 2)
		self.anim = animation.FuncAnimation(self.fig, self.animate, init_func=self.animate_init, frames=int(self.to/self.ts), interval=self.ts*1000, repeat=False, blit=True)
		self.writer = animation.AVConvWriter(fps=25)
		if save:
			self.anim.save(filename='tps.mp4', writer=self.writer)

	def animate_init(self):
		for p, pt, line,flabel, fv in zip(self.particles, self.pts, self.lines, self.force_labels, self.force_vectors):
			pt.set_data(p.pos[0], p.pos[1])
			pt.set_3d_properties(p.pos[2])
			line.set_data(p.pos[0], p.pos[1])
			line.set_3d_properties(p.pos[2])
			flabel.set_text("F=%.3fN\nlen=%.1fm" % (p.force, norm(p.pos)))
			flabel.set_position((p.pos[0] * 1.1, p.pos[1] * 1.1))
			flabel.set_3d_properties(p.pos[2] * 1.1)
			fv.set_3d_position([p.pos[0], p.Fdir[0]], [p.pos[1], p.Fdir[1]], [p.pos[2], p.Fdir[2]])
		return self.pts + self.lines + self.force_labels + self.force_vectors

	def animate(self, n):
		for p in self.particles:
			[p.add(po) for po in self.particles if po != p] # add other particle forces
		for p, pt, line, flabel, fv in zip(self.particles, self.pts, self.lines, self.force_labels, self.force_vectors):
			p.move(self.ts)
			pt.set_data(p.pos[0], p.pos[1])
			pt.set_3d_properties(p.pos[2])
			x = [0, p.pos[0]]
			y = [0, p.pos[1]]
			z = [0, p.pos[2]]
			line.set_data(x, y)
			line.set_3d_properties(z)
			flabel.set_text("F=%.3fN\nlen=%.1fm" % (p.force, norm(p.pos)))
			flabel.set_position((p.pos[0] * 1.1, p.pos[1] * 1.1))
			flabel.set_3d_properties(p.pos[2] * 1.1)
			fv.set_3d_position([p.pos[0], p.Fdir[0]], [p.pos[1], p.Fdir[1]], [p.pos[2], p.Fdir[2]])
		self.fig.canvas.draw()
		return self.pts + self.lines + self.force_labels + self.force_vectors

if __name__ == "__main__":
	sim = tpsim(6)
	show()