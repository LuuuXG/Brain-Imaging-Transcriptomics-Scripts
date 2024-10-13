require('crayon')

error <- red $ bold
warn <- yellow $ underline
note <- blue $ bold
dividing_line <- bgBlack$white$blurred('------------------------') %+%
  bgBlack$white('This is a dividing line') %+%
  bgBlack$white$blurred('------------------------\n')
