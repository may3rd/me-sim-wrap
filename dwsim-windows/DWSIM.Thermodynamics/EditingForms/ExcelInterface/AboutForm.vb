Imports Cudafy
Imports Cudafy.Host
Imports Yeppp

Public Class AboutForm

    Private Sub Button1_Click(sender As Object, e As EventArgs) Handles Button1.Click
        Me.Close()
    End Sub

    Private Sub AboutForm_Load(sender As Object, e As EventArgs) Handles MyBase.Load

        Application.EnableVisualStyles()

        Me.Text += " (" & My.Application.Info.DirectoryPath & "\DWSIM.xll)"

        Version.Text = "Excel Add-In Version " & My.Application.Info.Version.ToString

        Copyright.Text = My.Application.Info.Copyright

        LblOSInfo.Text = My.Computer.Info.OSFullName & ", Version " & My.Computer.Info.OSVersion & ", " & My.Computer.Info.OSPlatform & " Platform"
        LblCLRInfo.Text = "Microsoft .NET Framework, Runtime Version " & System.Runtime.InteropServices.RuntimeEnvironment.GetSystemVersion.ToString()
    
        Lblcpuinfo.Text = "Retrieving CPU info..."

        Lblcpusimd.Text = "Querying CPU SIMD capabilities..."

    End Sub

End Class