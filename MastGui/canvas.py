from Tkinter import *

class CanvasDemo(Frame):
    def __init__(self, width=200, height=200):
        Frame.__init__(self, root)
        self.canvas = Canvas(self)
        self.canvas.pack(fill="both", expand="1")
        self.canvas.create_rectangle(50, 25, 150, 75, fill="bisque", tags="r1")
        self.canvas.create_line(0,0, 50, 25, arrow="last", tags="to_r1")
        self.canvas.bind("<B1-Motion>", self.move_box)
        self.canvas.bind("<ButtonPress-1>", self.start_move)

    def move_box(self, event):
        deltax = event.x - self.x
        deltay = event.y - self.y
        self.canvas.move("r1", deltax, deltay)
        coords = self.canvas.coords("to_r1")
        coords[2] += deltax
        coords[3] += deltay
        self.canvas.coords("to_r1", *coords)
        self.x = event.x
        self.y = event.y

    def start_move(self, event):
        self.x = event.x
        self.y = event.y

root = Tk()
canvas = CanvasDemo(root)
canvas.pack()
mainloop()
