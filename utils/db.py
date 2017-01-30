import sqlite3
import os, sys
import datetime, time


class Tuning:
    def __init__(self, pars):
        self.time = pars[0]
        self.id = pars[1]
        self.wl = pars[2]
        self.charge = pars[3]
        self.comment = pars[4]

class ActionResult:
    def __init__(self, pars):
        self.tuning_id = pars[0]
        self.id = pars[1]
        self.sase_start = pars[2]
        self.sase_end = pars[3]

class ActionParameters:
    def __init__(self, pars):
        self.tuning_id = pars[0]
        self.action_id = pars[1]
        self.par_name = pars[2]
        self.start_value = pars[3]
        self.end_value = pars[4]


class PerfDB:
    def __init__(self, dbname="flash.db"):
        print("connecting to database ... ")
        self.db = sqlite3.connect(dbname)
        
    def new_tuning(self, params):
        print ('creating new tuning', params)
        with self.db:
            cursor = self.db.cursor()
            cursor.execute("insert into TUNINGS(TIME,CHARGE,WL, COMMENT) VALUES(?,?,?,?)",
                           (datetime.datetime.now(), params['charge'], params['wl'], "test"))
        
    def get_tunings(self):
        cursor = self.db.cursor()
        cursor.execute("select * from TUNINGS")
        return [Tuning(r) for r in cursor.fetchall()]

    def current_tuning_id(self):
        return self.get_tunings()[-1].id  # TODO: not effective

    
    def new_action(self, tuning_id, start_sase, end_sase):
        print ('creating new action')
        with self.db:
            cursor = self.db.cursor()
            cursor.execute("insert into ACTIONS(TUNING_ID,SASE_START,SASE_END) VALUES(?,?,?)",(tuning_id, start_sase,end_sase))


    def get_actions(self, tuning_id = None):
        cursor = self.db.cursor()
        if tuning_id == None: tuning_id = self.current_tuning_id()
        cursor.execute("select * from ACTIONS WHERE TUNING_ID=:Id", {'Id': tuning_id})
        return [ActionResult(r) for r in cursor.fetchall()]

    def current_action_id(self):
        return self.get_actions()[-1].id  # TODO: not effective


    def add_action_parameters(self, tuning_id, action_id, param_names, start_vals, end_vals):
        print('updating action', tuning_id, action_id)
        with self.db:
            cursor = self.db.cursor()

            for i in range(len(param_names)):
                cursor.execute("insert into PARAMETERS(TUNING_ID,ACTION_ID, PAR_NAME,PAR_START_VALUE,PAR_END_VALUE) VALUES(?,?,?,?,?)",
                               (tuning_id, action_id, param_names[i], start_vals[i], end_vals[i]))

    def get_action_parameters(self, tuning_id, action_id):
        cursor = self.db.cursor()
        cursor.execute("select * from PARAMETERS WHERE TUNING_ID=:tid and ACTION_ID = :aid", {'tid': tuning_id, 'aid': action_id})
        #return cursor.fetchall()
        return [ActionParameters(r) for r in cursor.fetchall()]

    def add_machine_parameters(self, tuning_id, params):
        print ('updating machine parameters for tuning ', tuning_id)
        with self.db:
            cursor = self.db.cursor()
            for k in params.keys():
                cursor.execute("insert into MACHINE_STATE(TUNING_ID,PAR_NAME,PAR_VALUE) VALUES(?,?,?)",
                               (tuning_id, k, params[k]))

    def get_machine_parameters(self, tuning_id):
        cursor = self.db.cursor()
        cursor.execute("select * from MACHINE_STATE WHERE TUNING_ID=:tid",{'tid':tuning_id})
        return cursor.fetchall()
        
    def close(self):
        self.db.close()

def create_db(dbname="flash.db"):

    db = sqlite3.connect(dbname)
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

def test_new_tunings(dbname):
    db = PerfDB(dbname=dbname)
    db.new_tuning({'wl':13.6, 'charge':0.1,'comment':'test tuning'}) # creates new tuning record (e.g. for each shift); 
    tunings = db.get_tunings()
    print ('current tunings', [(t.id, t.time, t.charge, t.wl) for t in tunings])
    tune_id = db.current_tuning_id()
    print ('current id', tune_id)
    db.add_machine_parameters(tune_id, params = {"hbar":1.0e-34, "nbunh":"ff"})
    print ('current machine parameters', db.get_machine_parameters(tune_id))


def test_new_action(dbname):
    db = PerfDB(dbname=dbname)
    tune_id = db.current_tuning_id()
    print ('new action for tune_id', tune_id)
    db.new_action(tune_id, start_sase = 1.0, end_sase = 150)
    print ('current actions in tuning', [(t.id, t.tuning_id, t.sase_start, t.sase_end) for t in db.get_actions()])

def test_add_action_parameters(dbname):
    db = PerfDB(dbname=dbname)
    tune_id = db.current_tuning_id()
    action_id = db.current_action_id()
    print ('updating', tune_id, action_id)
    db.add_action_parameters(tune_id, action_id, param_names = ["H1","Q1", "test"], start_vals = [0.1,0.2, "test"], end_vals=[1.1, 1.3, "gh"])
    print ('current actions', [(t.id, t.tuning_id, t.sase_start, t.sase_end) for t in db.get_actions()])
    print ('current action parameters', [(p.tuning_id, p.action_id, p.par_name, p.start_value, p.end_value) for p in db.get_action_parameters(tune_id, action_id)])



if __name__ == "__main__":
    # tests
    dbname = "test2.db"
    create_db(dbname)
    test_new_tunings(dbname)
    test_new_action(dbname)
    test_add_action_parameters(dbname)
    test_new_action(dbname)
    test_add_action_parameters(dbname)
