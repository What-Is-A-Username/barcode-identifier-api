from django.db.models import TextChoices

class DatabasePermissions(TextChoices): 
    '''
    Levels of permissions to access, edit and run a given database. These permissions hold regardless of whether the database is public or not.
        - CAN_EDIT_DB:      Allow user to edit database
        - CAN_RUN_DB:       Allow user to run BLAST on database 
        - CAN_VIEW_DB:      Allow user to view database data

    A superuser and the database owner can perform all actions for a given database.

    An unauthenticated user without any DatabasePermissions can only have view and run access, provided that the database is public. 

    Otherwise, A user marked with DENY_ACCESS will be denied access to all actions, regardless of the permissions specified above. 

    An authenticated user without any DatabasePermissions can only have view and run access, contingent that the database is public.

    If database is not specified (i.e. None), then permissions will return True if the action is permissible for any, but not necessary all, databases.
    
    '''

    # Explicitly prevent database from being viewed under any circumstances. Overrides 'public'
    DENY_ACCESS = 'deny_access'
    # Allow user to run BLAST on database 
    CAN_RUN_DB = 'can_run_db'
    # Allow user to view database. Overrides 'public'
    CAN_VIEW_DB = 'can_view_db'
    # Allow user to edit database data. Overrides 'public'
    CAN_EDIT_DB = 'can_edit_db'