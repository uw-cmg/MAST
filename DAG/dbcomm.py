import pymongo
from pymongo import MongoClient
import datetime
now = datetime.datetime.utcnow()

#def main():
#   connection = MongoClient()
connection = MongoClient('localhost',27017)
db =connection.test_database
collection = db.test_collection

post = {"author": "Mike",
        "text": "My first blog post!",
        "tags": [ "mongodb", "python", "pymongo"],
        "date": now()}
posts = db.posts
posts.insert(post)
db.collection_names()
collection = db.test_collection
posts.find_one()
print 'posts.find_one({"author":"Mike"})'
posts.find_one({"author":"Mike"})
print 'posts.find_one({"author":"Eliot"})'
posts.find_one({"author":"Eliot"})

post_id ='5143cbcf1310b3025487ffbe'
print "posts.find_one({'_id':post_id})"
posts.find_one({'_id':post_id})

# start mastdb in mangodb
mastdb = connection['mast-database']
connection = MongoClient('localhost',27017)
dishes = mastdb.dishes

dish={"dish_id":"d_1",
      "jobs_id":["j_1","j_2","j_3"],
      "submit":now(),
      "end_time":"NULL"}
dishes.insert(dish)

#find all jobs which are in a dish



            
