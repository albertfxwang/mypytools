#!/usr/bin/env python
import smtplib


def gmail(fromaddr, toaddrs, msg):
   username = 'astrokuang'
   password = 'YuHua1204'

   # The actual mail send
   server = smtplib.SMTP('smtp.gmail.com:587')
   server.starttls()
   server.login(username,password)
   server.sendmail(fromaddr, toaddrs, msg)
   server.quit()
   

