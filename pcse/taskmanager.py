# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), April 2014
"""Taskmanager which reads tasks from a database table and updates the status
when the tasks is finished.

Classes defined here:
 TaskManager
"""

import os
import logging
import sqlalchemy as sa
from sqlalchemy import select, and_, MetaData, Table
import socket

class TaskManager:
    """Class defines a taskmanager which reads from a table called 'tasklist'.

    Usage: tm = TaskManager(engine, dbtype=dbtype, tasklist)

    Public methods:
    get_task() - picks a 'Pending' task from the list
    set_task_finished(task) - set the task status to 'Finished'
    set_task_error(task) - set the task status to 'Error occurred'

    A task table could be created with the following SQL command
    (example taken from MySQL)::


        CREATE TABLE `tasklist` (
           `task_id` int(11) NOT NULL AUTO_INCREMENT,
           `status` char(16) DEFAULT NULL,
           `hostname` char(50) DEFAULT NULL,
           `process_id` int(11) DEFAULT NULL,
           `comment` varchar(200) DEFAULT NULL,
           `parameter1` int(11) DEFAULT NULL,
           `parameter2` decimal(10,2) DEFAULT NULL,

           ... Additional columns can be put here.

           PRIMARY KEY (`task_id`),
           KEY `status_ix` (`status`)
         );

    """
    validstatus = ['Pending', 'In progress', 'Finished',
                   'Error occurred']
    knowndatabases = ['mysql', 'oracle', 'sqlite']

#-------------------------------------------------------------------------------
    def __init__(self, engine, dbtype=None, tasklist='tasklist'):
        """Class constructor for TaskManager.
        
        Arguments:
        * engine - An SQLAlchemy engine object

        Keywords:
        * dbtype - the type of DB to connect either 'MySQL', 'ORACLE' or 'SQLite'
        * tasklist - Name of table to read tasks from, default 'tasklist'
        """
        db_ok = False
        if isinstance(dbtype, str):
            if dbtype.lower() in self.knowndatabases:
                db_ok = True
        if db_ok is False:
            msg = "keyword 'dbtype' should be one of %s" % self.knowndatabases
            raise RuntimeError(msg)
        if not isinstance(engine, sa.engine.base.Engine):
            msg = "Argument 'engine' should be SQLalchemy database engine, " \
                  "got %s" % engine
            raise RuntimeError(msg)

        self.dbtype = dbtype.lower()
        self.engine = engine
        self.logger = logging.getLogger("TaskManager")
        self.hostname = socket.gethostname()
        self.process_id = os.getpid()
        self.tasklist_tablename = tasklist

        # Check if tasklist exists and database is readable
        try:
            conn = self.engine.connect()
            metadata = MetaData(conn)
            self.table_tasklist = Table(tasklist, metadata, autoload=True)
        except Exception as e:
            msg = "Unable to connect or tasklist table doesn't exist!"
            self.logger.exception(msg)
            raise RuntimeError(msg)

#-------------------------------------------------------------------------------
    def get_task(self):
        "Return 'Pending' task to processing unit."
        
        conn = self.engine.connect()
        self._lock_table(conn)
        tasklist = self.table_tasklist
        s = select([tasklist], and_(tasklist.c.status=='Pending'), 
                   order_by=[tasklist.c.task_id],
                   limit=1)
        r = conn.execute(s)
        row = r.fetchone()
        r.close()
        if row is None:
            self._unlock_table(conn)
            return None
        else:
            task = dict(row)
            u = tasklist.update(tasklist.c.task_id==task["task_id"])
            conn.execute(u, status='In progress',  hostname=self.hostname,
                         process_id=self.process_id)
        self._unlock_table(conn)
        return task
    
#-------------------------------------------------------------------------------
    def _lock_table(self, connection):
        """Locks the TASKLIST table. Note that locking/unlocking is carried out
        by directly sending SQL to the database instead of using SQLAlchemy,
        which does not support this kind of table locking."""
    
        if self.dbtype=="mysql":
            connection.execute("LOCK TABLE %s WRITE" % self.tasklist_tablename)
        elif self.dbtype=="oracle":
            connection.execute("LOCK TABLE %s IN EXCLUSIVE MODE" %
                               self.tasklist_tablename)
        elif self.dbtype=="sqlite":
            pass # No locking needed for SQLite: assuming one client only.
        

#-------------------------------------------------------------------------------
    def _unlock_table(self, connection):
        "Unlocks the TASKLIST table"
    
        if self.dbtype=="mysql":
            connection.execute("UNLOCK TABLES")
        elif self.dbtype=="oracle":
            connection.execute("COMMIT")
        elif self.dbtype=="sqlite":
            pass # No locking needed for SQLite: assuming one client only.
        
#-------------------------------------------------------------------------------
    def set_task_finished(self, task, comment="OK"):
        "Sets a task to status 'Finished'"
        
        conn = self.engine.connect()
        self._lock_table(conn)
        u = self.table_tasklist.update(self.table_tasklist.c.task_id==task["task_id"])
        conn.execute(u, status='Finished', comment=comment)
        self._unlock_table(conn)
        
#-------------------------------------------------------------------------------
    def set_task_error(self, task, comment=None):
        "Sets a task to status 'Error occurred' with given comment"
        
        conn = self.engine.connect()
        self._lock_table(conn)
        u = self.table_tasklist.update(self.table_tasklist.c.task_id==task["task_id"])
        conn.execute(u, status='Error occurred', comment=comment)
        self._unlock_table(conn)
