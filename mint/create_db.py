import sqlite3
import os, sys

db = sqlite3.connect('flash.db')
cursor = db.cursor()

cursor.execute("drop table if exists TUNINGS")
sql = """CREATE TABLE TUNINGS (
         TIME  DATETIME NOT NULL,
         ID  INTEGER PRIMARY KEY,
         WL FLOAT,
         CHARGE FLOAT,
         COMMENT  CHAR(20))"""

cursor.execute(sql)

cursor.execute("drop table if exists ACTIONS")
sql = """CREATE TABLE ACTIONS (
         TUNING_ID INTEGER,
         ACTION_ID INTEGER PRIMARY KEY,
         SASE_START FLOAT,
         SASE_END FLOAT)"""

cursor.execute(sql)

cursor.execute("drop table if exists PARAMETERS")
sql = """CREATE TABLE PARAMETERS (
         TUNING_ID INTEGER,
         ACTION_ID INTEGER,
         PAR_NAME CHAR(20),
         PAR_START_VALUE FLOAT,
         PAR_END_VALUE FLOAT,
         PRIMARY KEY(TUNING_ID, ACTION_ID, PAR_NAME) )"""

cursor.execute(sql)

cursor.execute("drop table if exists MACHINE_STATE")
sql = """CREATE TABLE MACHINE_STATE (
         TUNING_ID INTEGER,
         PAR_NAME CHAR(20),
         PAR_VALUE FLOAT)"""

cursor.execute(sql)



db.close()