import threading
import time 
from queue import Queue
import numpy as np

print_lock = threading.Lock()

class MyThread(threading.Thread):
    def __init__(self, in_queue, out_Queue, args=(), kwargs=None):
        threading.Thread.__init__(self, args=(), kwargs=None)
        self.queue = in_queue
        self.daemon = True
        self.receive_messages = args[0]
        self.out_Queue = out_Queue

    def run(self):
        print(threading.currentThread().getName(), self.receive_messages)
        val = self.queue.get()
        val1 = self.queue.get()
        self.do_thing_with_message(val, val1)

    def do_thing_with_message(self, message, message1):
        if self.receive_messages:
            with print_lock:
                print (threading.currentThread().getName(), "Received")
                self.out_Queue.put(np.array(message) + np.array(message1))
                self.out_Queue.put(np.array(message1))

if __name__ == '__main__':
    threads = []
    outs = []
    for t in range(10):
        q_in = Queue()
        q_out = Queue()
        outs.append(q_out)
        threads.append(MyThread(q_in, q_out, args=(True,)))
        threads[t].start()
        time.sleep(0.1)

    for i, t in enumerate(threads):
        t.queue.put([10, 20])
        t.queue.put([5, 10])
        print(outs[i].get())
        print(outs[i].get())

    for t in threads:
        t.join()

'''
class test:
    num = 8
    bools = []
    def __init__(self):
        for i in range(self.num):
            print("starting thread", i)
            self.bools.append(0)
            t = threading.Thread(target=self.thread_func, args=(i,))
            t.daemon = True
            t.start()

        for i in range(len(self.bools)):
            self.bools[i] = 1

        time.sleep(2)

        for i in range(len(self.bools)):
            self.bools[i] = 2

    def thread_func(self, id):
        while True:
            if self.bools[id] == 1:
                print(f"Thread {id} here")
                self.bools[id] = 0
            
            if self.bools[id] == 2: return

inst = test()
'''