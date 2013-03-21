import pymongo
from pymongo import MongoClient
import datetime
def main():
%    connection = MongoClient()
    connection = MongoClient('localhost',27017)
    db =connection.test_database
    collection = db.test_collection

    post = {"author": "Mike",
            "text": "My first blog post!",
            "tags": [ "mongodb", "python", "pymongo"],
            "date": datetime.datetime.utcnow()}
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
    print 'posts.find_one({'_id':post_id})'
    posts.find_one({'_id':post_id})
