import os
import threading
import wx

from nrpslib import KnownCodes, get_adenylation_domains


wildcard = "Python source (*.py)|*.py|" \
            "All files (*.*)|*.*"


class LibThread(threading.Thread):
    def __init__(self, window, infile, known):
        super(LibThread, self).__init__()

        # Attributes
        self.window = window
        self.known = known
        self.infile = infile

    def run(self):
        known, domains = make_library(self.known, self.infile)
        oligos = ''
        for dom in domains:
            dom[0].make_library(known)
            dom[0].make_oligos()
            oligos = dom[0].write_oligo_string(oligos)
        wx.CallAfter(setattr, self.window, 'known', known)
        wx.CallAfter(setattr, self.window, 'domains', domains)
        wx.CallAfter(self.window.oligos.SetValue, oligos.strip()) #'\n'.join([str(dom[0]) for dom in domains]))
        wx.CallAfter(self.window._stop_busy)
        wx.CallAfter(self.window.saveFileDlgBtn.Enable)

def make_library(known, infile):
    """Calculate Fibonacci numbers
    using slow recursive method to demonstrate
    blocking the UI.
    """
    if not known:
        known = KnownCodes('lib/knowncodes.fasta')

    domains = get_adenylation_domains(infile, known)

    return known, domains


class MyApp(wx.App):
    def OnInit(self):
        self.frame = BoxSizerFrame(None, title="Nonribosomal wonder drugs")
        self.SetTopWindow(self.frame)
        self.frame.Show()

        return True


class BoxSizerFrame(wx.Frame):
    def __init__(self, *args, **kwargs):
        super(BoxSizerFrame, self).__init__(*args, **kwargs)

        # Attributes
        self.panel = BoxSizerPanel(self)

        # Layout
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.panel, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetInitialSize()


class BoxSizerPanel(wx.Panel):
    def __init__(self, *args, **kwargs):
        super(BoxSizerPanel, self).__init__(*args, **kwargs)

        # Attributes
        self.known = ''
        self.domains = ''
        self.current_directory = os.getcwd()
        self.infile = wx.TextCtrl(self)
        self.infile.Disable()

        self.openFileDlgBtn = wx.Button(self, label="Open fasta file")
        self.openFileDlgBtn.Bind(wx.EVT_BUTTON, self.on_open_file)

        self.checkbox = wx.CheckBox(self, label="Leading strand")
        self.checkbox.SetValue(wx.CHK_CHECKED)

        self.findDomainsBtn = wx.Button(self, label="Make oligos")
        self.findDomainsBtn.Bind(wx.EVT_BUTTON, self.find_domains)
        self.findDomainsBtn.Disable()
        self.oligos = wx.TextCtrl(self)

        self.saveFileDlgBtn = wx.Button(self, label="Save oligos")
        self.saveFileDlgBtn.Bind(wx.EVT_BUTTON, self.on_save_file)
        self.saveFileDlgBtn.Disable()

        self.timer = wx.Timer(self)
        self.prog = wx.Gauge(self)
        self.Bind(wx.EVT_TIMER, self._on_pulse, self.timer)


        # Layout
        self._DoLayout()

    def on_open_file(self, event):
        """Create and show the Open FileDialog"""
        dlg = wx.FileDialog(
            self, message="Choose a file",
            defaultDir=os.getcwd(),
            defaultFile="",
            wildcard=wildcard,
            style=wx.OPEN | wx.CHANGE_DIR
            )
        if dlg.ShowModal() == wx.ID_OK:
            self.infile.write(dlg.GetPath())
            self.findDomainsBtn.Enable()
        dlg.Destroy()

    def on_save_file(self, event):
        """
        Create and show the Save FileDialog
        """
        dlg = wx.FileDialog(
            self, message='Save file as ...',
            defaultDir=self.current_directory,
            defaultFile='out.txt', wildcard=wildcard, style=wx.SAVE
            )
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            with open (path, 'w') as f:
                f.write(self.oligos.GetValue())
        dlg.Destroy()

    def find_domains(self, event):
        """Find domains"""
        infile = self.infile.GetValue()
        #self.output.SetValue("") # clear output
        self._start_busy() # give busy feedback
        # Non-Blocking mode
        task = LibThread(self, infile, self.known)
        task.start()

    def _start_busy(self):
        self.timer.Start(100)
        self.findDomainsBtn.Disable()
        self.openFileDlgBtn.Disable()

    def _stop_busy(self):
        self.timer.Stop()
        self.prog.SetValue(0)
        self.findDomainsBtn.Enable()
        self.openFileDlgBtn.Enable()

    def _on_pulse(self, event):
        self.prog.Pulse() # Pulse busy feedback

    def _DoLayout(self):
        """Layout the controls"""
        vsizer = wx.BoxSizer(wx.VERTICAL)
        gsizer = wx.GridBagSizer(vgap=8, hgap=8)

        title = wx.StaticText(self, label="Nonribosomal oligo designer")
        infile_lbl = wx.StaticText(self, label="Sequence:")
        detail_lbl = wx.StaticText(self, label="Find domains:")

        # Title
        vsizer.Add(title, 0)
        vsizer.Add(self.prog, 0, wx.EXPAND)

        # Sequence file
        gsizer.Add(infile_lbl, (0, 1))
        gsizer.Add(self.openFileDlgBtn, (0, 2))
        gsizer.Add(self.infile, (0, 3), (1, 14), wx.EXPAND)
        gsizer.Add(self.checkbox, (0, 17))

        # Adenylation detection
        gsizer.Add(detail_lbl, (1, 1))
        gsizer.Add(self.findDomainsBtn, (1, 2))
        gsizer.Add(self.oligos, (2, 1), (10, 45), wx.EXPAND)

        # Sequence file
        #gsizer.Add(outfile_lbl, (0, 1))
        gsizer.Add(self.saveFileDlgBtn, (1, 3))

        # Add a spacer to pad out the right side
        # gsizer.Add((5, 5), (2, 17))
        # And another to the pad out the bottom
        # gsizer.Add((5, 5), (7, 0))
        vsizer.Add(gsizer, 0, wx.EXPAND|wx.ALL, 10)
        self.SetSizer(vsizer)



if __name__ == "__main__":
    app = MyApp(False)
    app.MainLoop()