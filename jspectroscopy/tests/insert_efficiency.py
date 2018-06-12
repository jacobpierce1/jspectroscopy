import jspectroscopy as spec
import time 

test_db = spec.spectrum_db( 'test_db', None, ( 32, 32 ), [ 'a', 'a' ] )

start_time = time.time() 

for i in range( 100 ) :

    test_db.insert_fit_data( -1, 31, 31, 0, None ) 

stop_time = time.time()

print( start_time - stop_time ) 
