#!/bin/sh

#  rectangle.sh
#
#
#  Created by Eric Tovar on 2/26/2018.
#
/Users/eric/Library/Mobile\ Documents/com~apple~CloudDocs/Texas\ A\&M/Math/spring17NumMthdsLab/triangle/triangle -npq28ea0.1 rectangle.poly
/Users/eric/Library/Mobile\ Documents/com~apple~CloudDocs/Texas\ A\&M/Math/spring17NumMthdsLab/triangle/triangle -rnpq28ea0.05  rectangle.1.poly
/Users/eric/Library/Mobile\ Documents/com~apple~CloudDocs/Texas\ A\&M/Math/spring17NumMthdsLab/triangle/triangle -rnpq28ea0.025  rectangle.2.poly
/Users/eric/Library/Mobile\ Documents/com~apple~CloudDocs/Texas\ A\&M/Math/spring17NumMthdsLab/triangle/triangle -rnpq28ea0.0125 rectangle.3.poly
/Users/eric/Library/Mobile\ Documents/com~apple~CloudDocs/Texas\ A\&M/Math/spring17NumMthdsLab/triangle/triangle -rnpq28ea0.00625 rectangle.4.poly
/Users/eric/Library/Mobile\ Documents/com~apple~CloudDocs/Texas\ A\&M/Math/spring17NumMthdsLab/triangle/triangle -rnpq28ea0.003125 rectangle.5.poly
/Users/eric/Library/Mobile\ Documents/com~apple~CloudDocs/Texas\ A\&M/Math/spring17NumMthdsLab/triangle/triangle -rnpq28ea0.0015625 rectangle.6.poly

# create FEM files here
triangle_2_fem rectangle.1
triangle_2_fem rectangle.2
triangle_2_fem rectangle.3
triangle_2_fem rectangle.4
triangle_2_fem rectangle.5
triangle_2_fem rectangle.6
triangle_2_fem rectangle.7
