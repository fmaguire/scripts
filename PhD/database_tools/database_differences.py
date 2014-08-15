#!/opt/anaconda/bin/python
import MySQLdb
import sys
import glob
import os

cwd = os.getcwd()

blastdb_folder = sys.argv[1]


def glob_and_trim(glob_string):
    untrimmed_list = glob.glob(glob_string)
    replace_string = glob_string.replace('*','')
    trimmed_set = {x.replace(replace_string, '') for x in untrimmed_list}
    return trimmed_set

def write_set(file_name, set_name):
    with open(file_name, 'wb') as fh:
         = list(set_name).sort()
        for item in set_name:
            fh.write(item + '\n')


os.chdir(blastdb_folder)
blastdb_fas = glob_and_trim("*.fas")
blastdb_phr = glob_and_trim("*.fas.phr")
blastdb_pin = glob_and_trim("*.fas.pin")
blastdb_psq = glob_and_trim("*.fas.psq")

blastdb_set = set.intersection(blastdb_fas,
                               blastdb_phr,
                               blastdb_pin,
                               blastdb_psq)

# present in blastdb folder but not in species
# list from MySQLdb


db_config = {"host": "",
             "db": "new_proteins",
             "user": "orchard",
             "passwd": ""}

con = MySQLdb.connect(host=db_config['host'],
                          user=db_config['user'],
                          passwd=db_config['passwd'],
                          db=db_config['db'])
cur = con.cursor()

mysql_cmd = "SELECT DISTINCT species FROM cider"
cur.execute(mysql_cmd)
mysql_species_list = cur.fetchall()
mysql_species_list_flat = [x[0] for x in mysql_species_list]
mysql_species_list_flat__ = {x.replace(' ', '_') for x in mysql_species_list_flat}


blastdb_not_mysql = blastdb_set.difference(mysql_species_list_flat__)

# present in MySQLdb but not in blastdb folder
mysql_not_blastdb = mysql_species_list_flat__.difference(blastdb_set)

os.chdir(cwd)
print 'in_blastdb_but_not_mysql'
print blastdb_not_mysql
print 'in_mysql_but_not_blastdb'
print mysql_not_blastdb


write_set('in_blastdb_but_not_mysql', blastdb_not_mysql)
write_set('in_mysql_but_not_blastdb', mysql_not_blastdb)
