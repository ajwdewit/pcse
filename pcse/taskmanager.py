"""Taskmanager which reads tasks from a database table and updates the status
when the tasks is finished.

Classes defined here:
 TaskManager
"""

from sqlalchemy import *

class TaskManager:
    """Class defines a taskmanager which reads from a table called 'tasklist'.

    Usage: tm = TaskManager(metadata, connection, dbtype=dbtype,
                            logger=logger, hostname=hostname)

    Public methods:
    get_task() - picks a 'Pending' task from the list
    set_task_finished(task) - set the task status to 'Finished'
    set_task_error(task) - set the task status to 'Error occurred'
    """
    
#-------------------------------------------------------------------------------
    def __init__(self, metadata, connection, dbtype=None, logger=None,
                 hostname=None, tasklist='tasklist'):
        """Class constructor for TaskManager.
        
        Arguments:
        * metadata - An SQLAlchemy database metadata object
        * connection - An SQLAlchemy database connection object
        
        Keywords:
        * dbtype - the type of DB to connect either 'MySQL', 'ORACLE', or 'SQLite'
        * hostname - The hostname of the machine executing the tasks
        * tasklist - Name of table to read tasks from, default 'tasklist'
        """
        
        self.dbtype = dbtype.lower()
        self.connection = connection
        self.logger = logger
        self.hostname = hostname
        self.tasklist_tablename = tasklist
        self.validstatus = ['Pending', 'In progress', 'Finished',
                            'Error occurred']
        self.knowndatabases = ['mysql', 'oracle','sqlite']
        
        # Check if database is known
        if self.dbtype not in self.knowndatabases:
            errstr =  "Unknown databases type!" +\
                      "Currently supported databases are: " +\
                      str(self.knowndatabases)
            raise RuntimeError(errstr)
        # Check if tasklist exists and database is readable
        try:
            self.table_tasklist = Table(tasklist, metadata, autoload=True)
        except Exception, e:
            raise e("Unable to connect or tasklist table doesn't exist!")
        
        
#-------------------------------------------------------------------------------
    def get_task(self):
        "Return 'Pending' task to processing unit."
        
        self._lock_table()
        tasklist = self.table_tasklist
        s = select([tasklist], and_(tasklist.c.status=='Pending'), 
                   order_by=[tasklist.c.task_id],
                   limit=1)
        r = self.connection.execute(s)
        row = r.fetchone()
        r.close()
        if row is None:
            self._unlock_table()
            return None
        else:
            task = dict(row)
            u = tasklist.update(tasklist.c.task_id==task["task_id"])
            self.connection.execute(u, status='In progress',
                                    hostname=self.hostname)
        self._unlock_table()
        return task
    
#-------------------------------------------------------------------------------
    def _lock_table(self):
        """Locks the TASKLIST table. Note that locking/unlocking is carried out
        by directly sending SQL to the database instead of using SQLAlchemy,
        which does not support this kind of table locking."""
    
        if self.dbtype=="mysql":
            self.connection.execute("LOCK TABLE %s WRITE" % self.tasklist_tablename) 
        elif self.dbtype=="oracle":
            self.connection.execute("LOCK TABLE %s IN EXCLUSIVE MODE" % \
                                    self.tasklist_tablename)
        elif self.dbtype=="sqlite":
            pass # No locking needed for SQLite: assuming one client only.
        

#-------------------------------------------------------------------------------
    def _unlock_table(self):
        "Unlocks the TASKLIST table"
    
        if self.dbtype=="mysql":
            self.connection.execute("UNLOCK TABLES")
        elif self.dbtype=="oracle":
            self.connection.execute("COMMIT")
        elif self.dbtype=="sqlite":
            pass # No locking needed for SQLite: assuming one client only.
        
#-------------------------------------------------------------------------------
    def set_task_finished(self, task):
        "Sets a task to status 'Finished'"
        
        self._lock_table()
        u = self.table_tasklist.update(self.table_tasklist.c.task_id==task["task_id"])
        self.connection.execute(u, status='Finished')
        self._unlock_table()
        
#-------------------------------------------------------------------------------
    def set_task_error(self, task):
        "Sets a task to status 'Error occurred'"
        
        self._lock_table()
        u = self.table_tasklist.update(self.table_tasklist.c.task_id==task["task_id"])
        self.connection.execute(u, status='Error occurred')
        self._unlock_table()
