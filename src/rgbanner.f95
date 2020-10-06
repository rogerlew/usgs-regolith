subroutine rgbanner(vrsn,bldate)
  character (len=7):: vrsn
  character (len=11):: bldate
  write (*,*) ''
  write (*,*) 'Regolith: Estimtate of regolith depth using'
  write (*,*) '   topographic data and geomorphic models'
  write (*,*) '       Version ', vrsn,', ',bldate
  write (*,*) '    By Rex L. Baum'
  write (*,*) '       U.S. Geological Survey'
  write (*,*) '-----------------------------------------------'
  write (*,*) 'Portions of this program include material'
  write (*,*) '     copyrighted (C) by'
  write (*,*) '      Absoft Corporation 1988-2017.'
  write (*,*) ''       
end subroutine rgbanner
