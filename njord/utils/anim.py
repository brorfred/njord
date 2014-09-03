import os
import tempfile

import pylab as pl
import figpref
class Movie:

    def __init__(self):
        self.reset()


    def reset(self):
        self.ffmpgF = ("-mbd rd -flags +mv4+aic -trellis 2 " +
                       " -cmp 2 -subcmp 2 -g 300 -y")
        self.oset = " -vcodec libx264 -qscale 1"
        self.tdir = tempfile.mkdtemp()
        print " === Images for video saved at " + self.tdir
        self.n = 0
        pl.close('all')
        figpref.current()

    def image(self,dpi=150):
        pl.savefig(self.tdir + '/mov_%06i.png' % self.n, dpi=dpi)
        self.n += 1

    def video(self, fn="movie.mp4",r=10):
        os.system('ffmpeg %s -r %i -i "%s/%s" %s %s' %
                  (self.ffmpgF, r, self.tdir,"mov_%06d.png", self.oset, fn) )
        print "---Movie saved as---"
        print "%s/%s" % (os.getcwd(), fn)

