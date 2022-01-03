import threading, time

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