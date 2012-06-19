#!/usr/bin/env python
'A generic CGI script'

##   Copyright (C) 2012 University of Texas at Austin
##  
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##  
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##  
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
import os, sys

import cgi
import cgitb; cgitb.enable()  # for troubleshooting

import shutil

def picture(directory,figure):
    png = figure+".png"
    fig = os.path.join("Fig",figure+".vpl")
    
    if os.chdir(directory):
        print "Content-type: text/html\n"
        print "<html><body>Wrong directory \"%s\".</body></html>" % directory
    elif os.system("source env.sh && scons " + fig) or \
            os.system("source env.sh && vpconvert pen=gd fat=3 serifs=n bgcolor=b %s %s" % (fig,png)):
        print "Content-type: text/html\n"
        print "<html><body>Madagascar failure.</body></html>"
    else:    
        print "Content-type: image/png\n"
        shutil.copyfileobj(open(png,"rb"), sys.stdout)

if __name__ == "__main__":
    form = cgi.FieldStorage()

    d = form.getvalue("dir")
    f = form.getvalue("fig")

    if not d or not f:
        print "Content-type: text/html\n"
        print "<html><body>Need dir= and fig=</body></html>"
        sys.exit(1)

    picture(d,f)
