function interrupt = keyboard_interrupt()

interrupt = 0;
if exist('testkeypress')==3 if testkeypress(' ')
  interrupt = 1;
end; end
if exist('Keytest')==8 if Keytest.test(' ')
  interrupt = 1;
end; end

